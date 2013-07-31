# mangle.py
"""Basic mangle routines for identifying which points are in which polygons.

For more on Mangle:
    http://space.mit.edu/~molly/mangle/

To build the included cython mangle_utils code, which will result in a ~4x
speed up for polyid commands:
    python setup.py build_ext --inplace

Example usage:
    import mangle
    polyfile = 'geometry-boss7.ply' #(also accepts fits files)
    mng = mangle.Mangle(polyfile)
    # ... read in a file with ra/dec coordinates as a numpy array
    polyids = mng.get_polyids(ra,dec)
    
    # The areas and weights of the found polygons can be gotten in two ways:
    mng.areas[polyids] == mng.get_areas(ra,dec)
    mng.weights[polyids] == mng.get_weights(ra,dec)
    
    # new mangle instance with just the first 10 polygons
    first_mng = mng[:10]

Requires numpy > 1.0, and pyfits > 3.0.4 for loading fits files.
"""
# Author: Martin White (UCB)
# Written 17-Mar-2010
# Modified: 17-Mar-2010 (Basic bug fixes)
# 22-Mar-2010 (Track pixel area)
# 20-Apr-2010 (Add setweight and writeply methods)
# 21-Apr-2010 (Add totalarea method)
# 12-May-2010 (Allow no space before "caps")
# 24-Jun-2010 jkp: (tweaks to gain ~40% speed improvement.)
# 29-Jun-2010 jkp: (numpy vectorized polyid for >10x speed improvement.)
# 22-Jul-2010 jkp: (completed and tested numpy version. Speed is ~1/2 that of the commandline mangle)
# 10-Aug-2010 jkp: (Added resorting of out-of-order polygons and altered original polylist structure)
# 10-Oct-2010 jkp: (Added __getitem__, returning a new Mangle instance given an index)
# 10-Mar-2011 jkp: (Added support for mangle_utils, with a 4x faster cythonized _incap)
# 18-May-2011 apw: (Added support using the photofield database tables in dr8db/spectradb)
# 10-Sep-2011 mecs: (allow ascii polygon files to be named .pol or .ply)
# 10-Sep-2011 mecs: (added option to keep existing id numbers from the polygon file)
# 15-Sep-2011 mecs: (added support to parse mangle header keywords in ascii and fits polygon files)
# 20-Sep-2011 mecs: (added method to write a polygon fits file)
# 27-Sep-2011 mecs: (added support to read in whatever extra column names are included in an input fits file)
# 27-Sep-2011 mecs: (added metadata structure to track of extra columns and methods to add and remove columns)
# 27-Sep-2011 mecs: (added support for reading and writing variable length array columns)
# 03-Nov-2011 mecs: (added support to read and write extra columns in ascii format)
# 29-Nov-2011 mecs: (fixed various bugs with importing ascii files with space before header keywords or negative weights/ids)
# 30-Nov-2011 mecs: (integrated with my graphmask plotting module)
# 21-Dec-2011 mecs: (added capability to use long doubles as in the real*10 version of mangle, added sortpolys function)
# 23-Dec-2011 mecs: (added which_pixel function and incorporated into get_polyids and polyid)
# 05-Jan-2011 jkp: Added write(), to intelligently write either .ply or .fits based on extension.
# 29-Feb-2012 ess: fixed bugs when mask not pixelized. Fixed bugs when 1-d arrays sent to get_polyids.
# 13-Aug-2012 mecs: Added append, rotate, and generate_shifted_copies functions
# 05-Jun-2013 jkp: removed implicit matplotlib dependency, by commenting out graphmask stuff.

import os
import re
import string
import copy
from itertools import izip
#import graphmask

import time

__author__ = "John Parejko, Martin White, Molly Swanson"
__version__ = "3.0 $Rev$"
__email__  = "john.parejko@yale.edu, mswanson@cfa.harvard.edu"

try:
    #raise ImportError # uncomment this to force pure-python code
    import mangle_utils
except ImportError as e:
    print "Error, could not import mangle_utils:",e
    print "Mangle routines will still function, but be slower."
    print "Check that you have built mangle_utils.so (see README.txt)."
    useUtils = False
else:
    useUtils = True

import numpy as np
# some of these will need to be called often, so explicitly name them
from numpy import array,sin,cos,fabs,pi,empty,zeros,ones,argsort,take,arange,ceil,floor
import pyfits

try:
    import longdouble_utils
except ImportError as e:
    uselongdoubles=False
    print "Error, could not import longdouble_utils:",e
    print "Mangle routines will use double precision for storing caps."
    print "To use longdouble precision (compatible with real*10 version of mangle), install mpmath:"
    print "http://code.google.com/p/mpmath/"
    print "and longdouble_utils:"
    print "[add link to longdouble_utils on github]"
    twopi = 2.*pi
else:
    uselongdoubles = True
    pi_long = np.arctan(np.longdouble(1))*4
    twopi = 2.*pi_long

class Mangle:
    """
    Implements basic mangle routines to identify polygons containing points.
    
    All functions assume ra,dec are in decimal degrees.

    Initialize with a single string containing the mask file, which can be
    either a polygon format file or a fits file.

    Default values, if the file does not specify them:
        weight == 0.
        area = -1.
        pixel = 0

    For working with large numbers of points, use the get_XXX() functions:
        get_polyids(), get_weight(), get_areas()
    which are up to 100x faster than the single value functions:
        polyid(), weight(), area()
    
    Supports slicing, e.g.:
        newmng = mng[:10]
        newmng2 = mng[mng.areas < 1e-7]
        newmng3 = mng[[1,3,5,7,9]]
    """

    def incap_spam(self,cap,x0,y0,z0):
        """
        Internal class function.

        Returns True if (theta,phi) lies in the cap specified by the
        4-vector "cap" containing (x,y,z,cm).
        """
        cd = 1.0-cap[0]*x0-cap[1]*y0-cap[2]*z0
        if cap[3]<0.0:
            if cd>fabs(cap[3]):
                return True
            else:
                return False
        else:
            if cd<cap[3]:
                return True
            else:
                return False

    def incap_vec(self,cap,x0,y0,z0):
        """
        Internal class function.

        Returns True for each (theta,phi) that lies in the cap specified by the
        4-vector "cap" containing (x,y,z,cm), and False for the rest."""
        cd = 1.0-cap[0]*x0-cap[1]*y0-cap[2]*z0
        return ((cap[3] < 0.0) & (cd>fabs(cap[3]))) | ((cap[3] > 0.0) & (cd<cap[3]))
    #...

    def inpoly_spam(self,polygon,theta,phi):
        """
        Internal class function.

        Returns True if (theta,phi) is in the polygon, i.e. if it is
        within all of the caps in the polygon.
        A polygon is a list giving the polygon number, the weight and
        then all of the caps (each cap is a 4-vector/4-list)."""
        # precompute the trig functions
        sintheta =  sin(theta)
        x0 = sintheta*cos(phi)
        y0 = sintheta*sin(phi)
        z0 = cos(theta)
        # this is faster than a for loop or list comprehension,
        # since it stops as soon as one returns false.
        return all(self.incap_spam(cap,x0,y0,z0) for cap in polygon)

    def inpoly_vec(self,polygon,x0,y0,z0):
        """
        Internal class function.

        Returns True for each (theta,phi) if it is in the polygon,
        i.e. if it is within all of the caps in the polygon.
        A polygon is a list giving the polygon number, the weight and
        then all of the caps (each cap is a 4-vector/4-list).
        """
        test = ones(len(x0),dtype=bool)
        if useUtils:
            if self.uselongdoubles:
                for cap in polygon:
                    #use long double version of _incap, called _incapl
                    test &= mangle_utils._incapl(cap,x0,y0,z0)
            else:
                for cap in polygon:
                    test &= mangle_utils._incap(cap,x0,y0,z0)
        else:
            for cap in polygon:
                test &= self.incap_vec(cap,x0,y0,z0)
        return test
    #...
    
    def which_simple_pixel(self,az,sinel):
        """Return the pixel ids for each pair of ra/dec for simple pixelization."""
        res = self.pixelization[0]
        pixel_start = (4**res-1)/3
        sinel = np.atleast_1d(sinel)
        az = np.atleast_1d(az)
        n = ceil((1-sinel)/2*(2**res))-1
        n[sinel==1] = 0
        m = floor((az%twopi)/twopi*(2**res))
        base_pix = 2**res*n+m
        pix = pixel_start+base_pix
        return pix.astype(int)
    #...
    
    def which_sdss_pixel(self,az,sinel):
        return None
    #...

    def which_pixel(self,az,sinel):
        """
        Return the mangle pixel numbers for each pair of ra,dec.
        pixelization info is stored in self.pixelization.
        self.pixelization[0] is the resolution.
        self.pixelization[1] is the pixelization scheme.
        !!! NOTE: only scheme 's' is currently implemented.
        """
        ### Math for this comes from which_pixel in the mangle2.2/src directory
        if self.pixelization == None:
            raise TypeError('No pixelization defined in this mangle instance.')
        if self.pixelization[0] < 0:
            raise TypeError('This python mangle wrapper does not yet support polyid adaptively pixelized files.')
        if self.pixelization[1] == 's':
            pix = self.which_simple_pixel(az,sinel)
        elif self.pixelization[1] == 'd':
            # TBD: still need to add support for sdsspix scheme, adaptive pixelization
            pix = self.which_sdss_pixel(az,sinel)
            raise TypeError('Only the "simple" pixelization scheme s is currently implemented.')
        else:
            raise TypeError('Only the "simple" pixelization scheme s is currently implemented.')
        return pix
    #...
    
    def inpoly_vec_pixels(self,radecpix,goodpolys,x0,y0,z0):
        """
        Internal class function for pixelized masks.
        
        radecpix is the result of calling self.which_pixel().
        Sets the polygon that xyz0[i] is in to goodpolys[i].
        """
        inpix = np.intersect1d(radecpix,self.pixels)
        for pix in inpix:
            radecs_in_pixel = (radecpix==pix)
            goodpolys_slice = goodpolys[radecs_in_pixel]
            for poly in self.pixel_dict[pix]:
                test = self.inpoly_vec(self.polylist[poly],x0[radecs_in_pixel],y0[radecs_in_pixel],z0[radecs_in_pixel])
                goodpolys_slice[test] = poly
            goodpolys[radecs_in_pixel] = goodpolys_slice
    #...

    def get_polyids(self,ra,dec):
        """
        Return the ID numbers of the polygons containing each RA/Dec pair.

        ra/dec should be numpy arrays in decimal degrees.

        !!! don't use on unbalkanized polygon files - this will only return 1 id per ra/dec.
        """
        if self.uselongdoubles:
            # force an upconversion to longdouble, in case ra/dec are float32 or float64
            theta = pi_long/180. * (90.0-np.array(dec,ndmin=1,dtype='f16'))
            phi = pi_long/180. * np.array(ra,ndmin=1,dtype='f16')
        else:
            # force an upconversion to double, in case ra/dec are float32
            theta = pi/180. * (90.0-np.array(dec,ndmin=1,dtype='f8'))
            phi = pi/180. * np.array(ra,ndmin=1,dtype='f8')
        sintheta = sin(theta)
        x0 = sintheta*cos(phi)
        y0 = sintheta*sin(phi)
        z0 = cos(theta)
        goodpolys = -ones(len(ra),dtype=int)
        # If we have a pixelized mask, we can reduce the number of polygons 
        # we have to check for each object.
        if self.npixels > 0:
            radecpix = self.which_pixel(phi,z0)
            self.inpoly_vec_pixels(radecpix,goodpolys,x0,y0,z0)
        else:
            for i,poly in izip(self.polyids,self.polylist):
                test = self.inpoly_vec(poly,x0,y0,z0)
                goodpolys[test] = i
        return goodpolys
    #...

    def contains(self,ra,dec):
        """
        Returns 1 if the point is found in the mask 0 if not.

        parameters
        ----------
        ra: scalar or array
        dec: scalar or array

        output
        ------
        a numpy array same length as ra,dec.  The value is 1 if it is in the
        mask (poly id found), 0 if not.  
        """
        ply=self.get_polyids(ra,dec)
        return np.where(ply >= 0, 1, 0)
    #...

    def get_areas(self,ra,dec):
        """
        Return the areas of the polygons containing each RA/Dec pair.

        ra/dec should be numpy arrays in decimal degrees.

        Return value is in steradians.

        Polygon areas are taken directly from the input file.
        """
        polyids=self.get_polyids(ra,dec)
        return self.areas[polyids]
    #...

    def get_weights(self,ra,dec):
        """
        Return the weights of the polygons containing each RA/dec pair.

        ra,dec should be numpy arrays in decimal degrees.

        The weight is taken directly from the input file."""
        polyids=self.get_polyids(ra,dec)
        return self.weights[polyids]
    #...

    def polyid(self,ra,dec):
        """
        Return the polyid of the polygon containing (ra,dec).
        
        ra,dec should be in decimal degrees.
        """
        if self.uselongdoubles:
            theta = pi_long/180. * (90.0-dec)
            phi   = pi_long/180. * ra
        else:
            theta = pi/180. * (90.0-dec)
            phi   = pi/180. * ra
        ipoly = -1

        if self.npixels > 0:
            pix=self.which_pixel(phi,cos(theta))
            polys_in_pixel=(self.pixels==pix)
            for i,poly in izip(self.polyids[polys_in_pixel],self.polylist[polys_in_pixel]):
                if self.inpoly_spam(poly,theta,phi):
                    ipoly = i
                    break
        else:
            for i,poly in zip(self.polyids,self.polylist):
                if self.inpoly_spam(poly,theta,phi):
                    ipoly = i
                    break
        return(ipoly)

    def weight(self,ra,dec):
        """
        Return the weight of the polygon containing (ra,dec).
        
        ra,dec should be in decimal degrees.
        """
        if self.uselongdoubles:
            theta = pi_long/180. * (90.0-dec)
            phi   = pi_long/180. * ra
        else:
            theta = pi/180. * (90.0-dec)
            phi   = pi/180. * ra
        weight= -1
        for w,poly in zip(self.weights,self.polylist):
            if self.inpoly_spam(poly,theta,phi):
                weight = w
                break
        return(weight)

    def area(self,ra,dec):
        """
        Return the area of the polygon containing (ra,dec).
        
        ra,dec should be in decimal degrees.

        Return value is in steradians.
        """
        if self.uselongdoubles:
            theta = pi_long/180. * (90.0-dec)
            phi   = pi_long/180. * ra
        else:
            theta = pi/180. * (90.0-dec)
            phi   = pi/180. * ra
        area  = -1.0
        for a,poly in izip(self.areas,self.polylist):
            if self.inpoly_spam(poly,theta,phi):
                area = a
                break
        return(area)

    def totalarea(self):
        """
        Return the total area and the 'effective' area in the mask.

        total area is the sum of the areas of each polygon.
        effective area is the area weighted by the completeness.
        """
        tot_area = self.areas.sum()
        eff_area = (self.weights*self.areas).sum()
        return((tot_area,eff_area))

    def compute_sector_areas(self):
        """
        Compute the area of each sector, and store them as self.sector_areas
        
        Requires that the sector numbers exist in sector.
        """
        sector_areas = {}
        for a,s in izip(self.areas,self.sector):
            if s not in sector_areas:
                sector_areas[s] = 0
            sector_areas[s] += a
        self.sector_areas = np.zeros(len(self.areas),dtype=self.areas.dtype)
        for i in range(len(self.areas)):
            self.sector_areas[i] = sector_areas[self.sector[i]]
    #...

    def set_weight(self,polyid,weight):
        """Set the weight of polygon 'polyid' to weight.
        
        polyid can be an int, index array or boolean index array.
        """
        self.weights[polyid] = weight

    def set_area(self,polyid,area):
        """Set the area of polygon 'polyid' to area.
        
        polyid can be an int, index array or boolean index array.
        """
        self.area[polyid] = area

    def set_pixel(self,polyid,pixel):
        """Set the pixel ID of polygon 'polyid' to pixel.
        
        polyid can be an int, index array or boolean index array.
        """
        self.pixel[polyid] = pixel

    def set_allarea(self,area):
        """Set the area of all polygons to area."""
        self.areas.flat = area

    def set_allweight(self,weight):
        """Set the weight of all polygons to weight."""
        self.weights.flat = weight

    def set_allpixel(self,pixel):
        """Set the weight of all polygons to weight."""
        self.pixels.flat = pixel
        
    def point_to_radec(self,ra,dec,ra_start=0,dec_start=0):
        """
        take a pattern of polygons centered at ra,dec = ra_start,dec_start and shift them by given ra,dec
        """
        psi=(90.0+ra_start)*pi/180.0
        theta=(dec)*pi/180.0
        phi=(-90.0-(ra+ra_start))*pi/180.0
        self.rotatepolys(psi,theta,phi)
        
    def rotate(self,azn=192.859508, eln=27.128336, azp=122.932):
        """
        perform a coordinate transform rotation from old frame to new frame
        uses same azn,eln,azp notation as mangle command rotate
        azp=azimuth of old pole wrt new frame in degrees
        azn=azimuth of new pole wrt old frame in degrees
        eln=elevation of new pole wrt old frame in degrees
        new north pole: (az_new,el_new)=(0,90), (az_old,el_old)=(azn,eln)
        old north pole: (az_old_el_old)=(0,90), (az_new,el_new)=(azp,eln)
        
        Example: old=equatorial (J2000), new=galactic
        ra,dec of north galactic pole = 192.859508, 27.128336
        galactic longitude of north celestial pole = 122.932
        so azn=192.859508, eln=27.128336, azp=122.932        
        
        """
        psi=(90.0+azn)*pi/180.0
        theta=(90.0-eln)*pi/180.0
        phi=(90.0-azp)*pi/180.0
        rotmat=self.rotatepolys(psi,theta,phi)
        return rotmat
        
    def rotatepolys(self,psi,theta,phi):
        """
        rotate polys using Euler angles phi,theta,psi in radians
        """

        rotate1= np.array([ 
                 [cos(psi), sin(psi), 0.0],
                 [-sin(psi), cos(psi), 0.0 ],
                 [0.0, 0.0, 1.0]
                 ])
        rotate2= np.array([ 
                 [1.0, 0.0, 0.0],
                 [0.0, cos(theta), sin(theta)],
                 [0.0, -sin(theta), cos(theta)]
                 ])
        rotate3= np.array([ 
                 [cos(phi), sin(phi), 0.0],
                 [-sin(phi), cos(phi), 0.0 ],
                 [0.0, 0.0, 1.0]
                 ])
        rotmat=np.dot(rotate3,np.dot(rotate2,rotate1))

        #list comprehension version
        self.polylist= array([None]+ #combine output of list comprehension with "None" in order to get the correct array type
            [np.vstack( #stack together rotated xyz vectors and original cm values
                    [np.dot(rotmat,xyzcm[:,:-1].T), #extract xyz vectors of caps and rotate them
                    xyzcm[:,-1]]  #extract cm values
                    ).T   
                for xyzcm in self.polylist])[1:] #loop through all polygon [x y z cm] arrays in polylist

        ## one-loop version
        #for i in range(self.npoly):
         #   self.polylist[i][:,:-1]=np.dot(rotmat,self.polylist[i][:,:-1].T).T
        
        ## two-loop version
        #for i in range(self.npoly):    
            #for cap in range(0,self.ncaps[i]):
            #    xyz=self.polylist[i][cap][:-1]
            #    xyz=np.asarray((rotate2*rotate1*np.matrix(xyz).T).T)[0]
            #    self.polylist[i][cap][:-1]=xyz
        return rotmat

    def generate_shifted_copies(self,ra,dec,weights=None):
        """
        takes a mangle instance as a base element and returns a new 
        mangle instance of many copies of that element, rotated to 
        the center points specified in ra,dec. Element is assumed
        to be centered on ra,dec=0,0

        Examples: 
        Generate camera footprint with polygons for 1 ccd as
        the original element and ra,dec as a list of the angular locations
        of the ccds within the camera.
        camera_polys=ccd_poly.generate_shifted_copies(ccdpos_ra,ccdpos_dec)

        Generate a series of camera pointing footprints with the camera
        polygons as the polygon_element and ra,dec as the ra/dec
        centers of the target pointings.
        footprint_polys=camera_polys.generate_shifted_copies(target_ra,target_dec)

        ra,dec should be numpy arrays in decimal degrees
        """
        #TODO
        #take an optional recarray input of the same length as ra, dec and add the values
        #as additional columns in the mangle instance
        
        #function to make a copy of the element, rotate it, and return it
        def rotated_element(element,az,el,weight=None):
            rotated_element=copy.deepcopy(element)        
            rotated_element.point_to_radec(az,el)
            if weight is not None:
                rotated_element.set_allweight(weight)
            return rotated_element
        if weights is None:
            weights=np.repeat(None,len(ra))
        #initialize output polygons with the first rotated element
        elements=rotated_element(self,ra[0],dec[0],weights[0])       
        #iterate through ra,dec positions, and append on a new rotated element for each one
        [elements.append(element) for element in  (rotated_element(self,r,d,w) for r,d,w in zip(ra[1:],dec[1:],weights[1:])) ]
        return elements

    def write(self,filename,**kwargs):
        """
        Write these polygons to either a .fits or .ply/.pol file, depending on the 
        extension of filename.
        
        Passes any keyword arguments on to either writeply() or write_fits_file().
        """
        if '.fits' in filename:
            self.write_fits_file(filename,**kwargs)
        if '.ply' in filename:
            self.writeply(filename,**kwargs)
        if '.pol' in filename:
            self.writeply(filename,**kwargs)
    #...

    def writeply(self,filename,keep_ids=False,write_extra_columns=False,weight=None):
        """
        Write a Mangle-formatted polygon file containing these polygons.

        If keep_ids=True, polygon ids from self.ids will be written as the polygon ids.

        If a name in self.names is provided as the weight argument, output polygons
        will have their weights given by the values in that column.
        
        If write_extra_columns=True, supplemental files containing additional information
        will be written in this format:
        for a polygon file named poly.pol or poly_obstime.pol
        (to denote that the polygons are weighted by observation time)
        Additional columns will be written as poly.maglims, poly.obstime, poly.XXX where XXX
        is the name of the column.
        If a name is provided as the weight argument, e.g., weight='obstime', then that
        column will be written as the weights, and if 'obstime' is contained in the filename,
        e.g., polys_obstime.pol, additional columns will be written as poly.XXX (not poly_obstime.XXX).
        Additional column files will have number of rows equal to the number of polygons. 
        This allows alternate weights and extra information about the polygons that are
        stored in a fits file to be written in ascii format.        
        """
        ff = open(filename,"w")
        ff.write("%d polygons\n"%len(self.polylist))
        ff.write("real %d\n"%self.real)
        if self.pixelization is not None:
            ff.write("pixelization %d%s\n"%(self.pixelization[0],self.pixelization[1]))
        if self.snapped == True:
            ff.write("snapped\n")
        if self.balkanized == True:
            ff.write("balkanized\n")
        if ('ids' in self.names) & keep_ids==True:
            ids_to_write=self.ids
        else:
            ids_to_write=self.polyids
        for i in range(self.npoly):
            str = "polygon %10d ( %d caps,"%(ids_to_write[i],len(self.polylist[i]))
            if self.uselongdoubles:
                weightstr=longdouble_utils.longdouble2string(self.weights[i],n=19)
                areastr=longdouble_utils.longdouble2string(self.areas[i],n=19)
                str+= " %s weight, %d pixel, %s str):\n"%(weightstr,self.pixels[i],areastr)
            else:
                str+= " %.15g weight, %d pixel, %.15g str):\n"%(self.weights[i],self.pixels[i],self.areas[i])
            ff.write(str)
            if self.uselongdoubles:
                for cap in self.polylist[i]:
                    capstr=[longdouble_utils.longdouble2string(c,n=19) for c in cap]
                    ff.write(" %s %s %s %s\n"%\
                             (capstr[0],capstr[1],capstr[2],capstr[3]))                    
            else:
                for cap in self.polylist[i]:
                    ff.write(" %.15g %.15g %.15g %.15g\n"%\
                             (cap[0],cap[1],cap[2],cap[3]))

        ff.close()

        #write extra columns stored in this format:
        # polygon file named poly.pol or poly_obstime.pol (to denote that the polygons are weighted by observation time)
        # additional columns will be written as poly.maglims, poly.obstime, poly.XXX where XXX is the name of the column
        # if a name is provided as the weight argument, e.g., weight='obstime', then that column will be written as the weights, and if 'obstime'
        # is contained in the filename, e.g., polys_obstime.pol, additional columns will be written as poly.XXX (not poly_obstime.XXX)
        # additional column files will have number of rows equal to the number of polygons 
        # this allows alternate weights and extra information about the polygons that are stored in a fits file to be written in ascii format        
        if write_extra_columns:
            base_filename=string.split(filename,'.')[0]
            # if column to weight by is provided, e.g., weights='obstime', check to see if base filename ends
            # in _obstime and strip it off if it does before generating extra column filenames 
            if weight is not None:
                if weight==string.split(base_filename,'_')[-1]:
                    base_filename='_'.join(string.split(base_filename,'_')[:-1])
            for name in self.names:
                asciiformat=self.metadata[name]['asciiformat']
                if asciiformat is not None:
                    np.savetxt(base_filename+'.'+name, vars(self)[name],fmt=asciiformat)

    def __getitem__(self,idx):
        """Slice this mangle instance by the index array idx.
        
        Return a mangle instance containing only those polygons
        referenced in idx, which can be a single integer, array or list
        of indicies, a boolean index array, or a slice.
        """
        mng2=copy.deepcopy(self)
        mng2.polyids = self.polyids[idx]
        # slices also have no len(), so we need to test after indexing
        # incase a slice was passed.
        try:
            mng2.npoly = len(mng2.polyids)
        except TypeError:
            idx = [idx,]
            mng2.polyids = self.polyids[idx]
            mng2.npoly = len(mng2.polyids)
        mng2.polylist = mng2.polylist[idx]
        mng2.areas = mng2.areas[idx]
        mng2.weights = mng2.weights[idx]
        mng2.ncaps = mng2.ncaps[idx]
        mng2.pixels = mng2.pixels[idx]
        for name in mng2.names:
            vars(mng2)[name]=vars(mng2)[name][idx]
        return mng2
    #...

    def append(self,mng):
        """
        Append the polygons in mangle instance mng onto mangle instance self.
        Example: 
        polys1=mangle.Mangle('test1.pol')
        polys2=mangle.Mangle('test2.pol')
        polys1.append(polys2)       
        """
        assert self.pixelization == mng.pixelization, "pixelizations of mangle instances don't match"
        if self.pixelization is not None:
            assert self.pixelization[0] != -1, "can't append to a set of polygons pixelized adaptively"
        assert self.balkanized == mng.balkanized, "one set of polygons is balkanized and one is not"
        assert self.snapped == mng.snapped, "one set of polygons is snapped and one is not"
        assert self.real == mng.real, "floating point precision of mangle instances don't match"
        assert self.names == mng.names, "column names of mangle instances don't match"
        self.npoly += mng.npoly
        self.polyids = np.append(self.polyids,mng.polyids)
        self.areas = np.append(self.areas,mng.areas)
        self.weights = np.append(self.weights,mng.weights)
        self.ncaps = np.append(self.ncaps,mng.ncaps)
        self.pixels = np.append(self.pixels,mng.pixels)
        self.polylist = np.append(self.polylist,mng.polylist)
        for name in self.names:
            vars(self)[name]=np.append(vars(self)[name],vars(mng)[name])
    #...

    def sortpolys(self,sortbycolumn='polyids',kind='quicksort'):
        """
    Function to sort polygons in-place.
    Parameters:
    sortbycolumn: indicates the column of the polygon structure to use as
                  the sorting index. Default is to sort by polyids.
    kind: sorting algorithm to use.  Default is quicksort.
        """
        sorter = argsort(vars(self)[sortbycolumn],kind=kind)
        self.polylist = take(self.polylist,sorter)
        self.weights = self.weights[sorter]
        self.areas = self.areas[sorter]
        self.polyids = self.polyids[sorter]
        self.ncaps = self.ncaps[sorter]
        self.pixels = self.pixels[sorter]
        for name in self.names:
            vars(self)[name]=vars(self)[name][sorter]

    # def drawpolys(self,skymap=None,**kwargs):
    #     """
    # Function to draw a set of mangle polygons.  Wrapper for graphmask.plot_mangle_map
           
    # Return value(s): - matplotlib.collections.PolyCollection object containing
    #                    the polygons that were plotted
    #                  - (optional) Skymap instance defining the map projection
    #                  e.g. : p=plot_mangle_map('mypolys.pol')
    #                         # p is the plotted PolyCollection
    #                         p,m=plot_mangle_map('mypolys.pol')
    #                         # p is the plotted PolyCollection, m is the Skymap instance

    # Optional keyword arguments:

    # skymap: a Skymap instance defining the map projection to plot in.  Use this
    #         to plot additional polygons in the same plot, e.g.,
    #         polys1=mangle.Mangle('polys1.pol')
    #         polys2=mangle.Mangle('polys2.pol')
    #         p1,m=polys1.drawpolys()
    #         p2=polys2.drawpolys(skymap=m)

    # outfilename: name of image file to save.  if None, just draw on screen.
    #              Accepts any format valid for matplotlib's savefig
    #              (usually png, pdf, ps, eps, and svg)

    # graphicsfilename: name of graphics (.list format) file to save to speed
    #                   up subsequent plotting (with graphmask.plot_mangle_map).
    #                   Should have a .list suffix.
    #                   e.g. : mypolys=mangle.Mangle('mypolys.pol')
    #                          p,m=mypolys.drawpolys(graphicsfilename='mypolys.list')
    #                          #this takes a long time, and after looking at the plot you decide
    #                          #it would be much better to plot the polygons in red.
    #                          p=plot_mangle_map('mypolys.list',facecolor='r')
    #                          #this plots much faster.

    #                    Note that the .list file saved is projection-specific
    #                    (and plot-range-specific) since the polygons are trimmed
    #                    to the map projection region plotted.  If you change the
    #                    projection or the plot range, you might need to use the
    #                    original polygon file again.                      

    # pointsper2pi: number of points along each circle edge used in making the
    #               .list file with mangle's poly2poly.  default is 30.

    # cmap: colormap used to plot the weight values of the polygons.  Can be any
    #       valid matplotlib colormap, or a string defining one, e.g. cmap=cm.jet
    #       or cmap='jet'.  Default is 'gray_r' which plots black/grayscale polygons
    #       on a white background.

    # plottitle: a title for the plot.

    # autoscale: if autoscale is True, the plot range will be determined by the
    #            min and max azimuth and elevation of the polygons. Default is False.

    # enlarge_border: fraction to enlarge the range determined by autoscale to give
    #                 a reasonable border around the plot.  Default is .1 (10%)

    # bgcolor: background color of the map.  If 'auto', use the color at the zero
    #          value of the colormap.  Default is 'auto'.

    # drawgrid: if True, draw grid lines for azimuth and elevation (latitude and longitude)
    #           and labels for them.

    # gridlinewidth: line width for the grid lines.  Use gridlinewidth=0 to draw
    #                axis labels but no grid lines.

    # minaz,maxaz,minel,maxel: if these are defined, use these values to define the
    #                          corner regions of the plot.

    # cenaz: central azimuth value to use for map projection (equiv to Basemap's lon_0)

    # cenel: central elevation to use for map projeciton (equiv to Basemap's lat_0)

    # Additional keyword arguments are passed to the Skymap class, which accepts
    # any keyword arguments recognized by the Basemap class.  Commonly used ones are:

    # projection: the map projection used by Basemap.  print the basemap module variable
    #             ``supported_projections`` to see a list of the options.

    # width, height: the width and height of the region to be plotted,
    #                in units of degrees on the sky.  cenaz,cenel,width,height
    #                can be used instead of minaz,maxaz,minel,maxel to define
    #                the plotted region.

    # Further additional keyword arguments not accepted by Skymap/Basemap are passed
    # to the draw_weighted_polygons function, which accepts

    # holes: can be 'cw'==clockwise,'ccw'==counterclockwise, or None.  If 'cw' ('ccw'),
    #        polygons with the edge points ordered clockwise (counterclockwise) are
    #        treated as holes and plotted in the background color.  Mangle polygons
    #        that have holes in them will be written as separate polygons in the .list
    #        format, with the exterior border in counterclockwise order and the holes
    #        in clockwise order, so 'cw' is the default.
    #        If you know your polygons have no holes, using None can speed up
    #        the plotting.
    #        Note that clockwise and counterclockwise are defined as when looking
    #        up at the sky, so it is reversed from what geographers would expect
    #        from looking down at the surface of the Earth.

    # emptyweight: value of polygon weight indicating a polygon that should not be plotted.
    #              Default is 0, but emptyweight can be set to something else to allow for
    #              weights whose allowed value range includes 0.

    # draw_colorbar: if True, draw a colorbar showing the color scale of the weight
    #                values plotted.  Default is True unless all weight values are
    #                equal. If you want to control how the colorbar is drawn,
    #                use draw_colorbar=False and then use the returned
    #                PolyCollection object as your ScalarMappable object.

    # vmin,vmax: minimum and maximum values for the colormap scale.

    # Still further keywords are passed on to the PolyCollection plotter, which can
    # used to customize the appearance of the plotted polygons, e.g.,
    # facecolor='red',edgecolor='black',linewidth=0.5.

    # Examples:
    #  #plot overlapping SDSS plates in a small example patch of sky,
    #  #using a transparent face color to illustrate the overlap.
    #  dr9plates=mangle.Mangle('dr9plates.pol')
    #  dr9plates.drawpolys(minaz=117.5,maxaz=125.5,minel=30, maxel=36.5,facecolor=[0, 0, 1, .2],linewidth=0)
    #  #plot a boss geometry file starting from a .fits file, and save a graphics .list file for faster plotting later
    #  geometry=mangle.Mangle('boss_geometry_2011_05_20.fits')
    #  p,m=geometry.drawpolys(graphicsfilename='boss_geometry_2011_05_20.list',cenaz=270,gridlinewidth=.1,projection='moll')
    #  """
         
    #     if 'graphicsfilename' in kwargs:
    #         graphicsfilename=kwargs.pop('graphicsfilename')
    #     else:
    #         graphicsfilename=None
    #     if 'pointsper2pi' in kwargs:
    #         pointsper2pi=kwargs.pop('pointsper2pi')
    #     else:
    #         pointsper2pi=30
    #     if skymap is None:
    #         if hasattr(self,'graphicsfilename'):
    #             p,m=graphmask.plot_mangle_map(self.graphicsfilename,graphicsfilename=graphicsfilename,pointsper2pi=pointsper2pi,**kwargs)
    #         else:
    #             p,m=graphmask.plot_mangle_map(self,graphicsfilename=graphicsfilename,pointsper2pi=pointsper2pi,**kwargs)
    #     else:
    #         m=skymap
    #         if hasattr(self,'graphicsfilename'):
    #             azel,weight=m.get_graphics_polygons(self.graphicsfilename,graphicsfilename=graphicsfilename,pointsper2pi=pointsper2pi)
    #         else:
    #             azel,weight=m.get_graphics_polygons(self,graphicsfilename=graphicsfilename,pointsper2pi=pointsper2pi) 
    #         p=graphmask.draw_weighted_polygons(m.world2pix(azel),weight=weight,**kwargs)
    #     if graphicsfilename is not None:
    #         self.graphicsfilename=graphicsfilename
    #     howmany=graphmask.expecting()
    #     if howmany<2:
    #         return p
    #     else:
    #         return p,m
        
    ## def graphics(self):
    ##     """Return an array of edge points, to plot these polygons.

    ##     Calls the command-line 'poly2poly' to generate a temporary graphics
    ##     file and return the output.

    ##     The returned array has len() == npoly, with each element
    ##     being an array of ra/dec pairs. Plot a polygon with, e.g.:
    ##         polys = mng.graphics()
    ##         pylab.plot(polys[0]['ra'],polys[0]['dec'])
    ##     or all of them (this may take a while,if there are a lot of polygons):
    ##         polys = mng.graphics()
    ##         for poly in polys:
    ##             pylab.plot(poly['ra'],poly['dec'])
    ##     """
    ##     import tempfile
    ##     import subprocess
    ##     tempIn = tempfile.NamedTemporaryFile('rw+b')
    ##     self.writeply(tempIn.name)
    ##     tempOut = tempfile.NamedTemporaryFile('rw+b')
    ##     call = ' '.join(('poly2poly -q -ol30',tempIn.name,tempOut.name))
    ##     subprocess.call(call,shell=True)
    ##     # NOTE: poly2poly -ol always outputs a separate weight file.
    ##     tempIn.close()
    ##     # the mangle list format has polygons delimited by NaN pairs.
    ##     dtype=[('ra','>f8'),('dec','>f8')]
    ##     data=np.loadtxt(tempOut,dtype=dtype)
    ##     weights=np.loadtxt(tempOut.name+'.weight',dtype=dtype)
    ##     tempOut.close()
    ##     os.remove(tempOut.name+'.weight')
    ##     i = 0
    ##     polys = np.empty(len(weights),dtype='object')
    ##     temp = []
    ##     for x in data:
    ##         if np.isnan(x[0]) and np.isnan(x[1]):
    ##             polys[i] = np.array(temp,dtype=dtype)
    ##             temp = []
    ##             i += 1
    ##         else:
    ##             temp.append(x)
    ##     return polys
    #...

    def read_ply_file(self,filename,read_extra_columns=False,real=None):
        """
        Read in polygons from a .ply file.

        If read_extra_columns=True, supplemental files containing additional information
        will be read in this format:
        For a polygon file named poly.pol or poly_obstime.pol (to denote that the polygons
        are weighted by observation time), additional columns can be provided as poly.maglims,
        poly.obstime, poly.XXX where XXX is whatever tag applies to the column.
        An additional column will be added to the mangle polygon object named XXX. 
        Additional column files should have number of rows equal to the number of
        polygons in poly.pol. This allows alternate weights and extra information about the
        polygons that are stored in a fits file to be written in ascii format

        'real' keyword indicates whether to use doubles (real=8) or long doubles (real=10)
        as corresponding to the real*8 and real*10 versions of mangle. If real is None, it will
        default to real=8 unless specified otherwise in the polygon file.
        """
        import glob
        # It's useful to pre-compile a regular expression for a mangle line
        # defining a polygon.
        rePoly = re.compile(r"polygon\s+([-+]?\d+)\s+\(\s*([-+]?\d+)\s+caps")
        reWeight=re.compile(r"([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)\s+weight")
        reArea=re.compile(r"([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)\s+str")
        rePixel = re.compile(r"([-+]?\d+)\s+pixel")
        rePixelization = re.compile(r"pixelization\s+([-+]?\d+)([sd])")
        rePolycount = re.compile(r"([-+]?\d+)\s+polygons")
        reSnapped = re.compile(r"snapped")
        reBalkanized = re.compile(r"balkanized")
        reReal = re.compile(r"real\s+([-+]?\d+)")
        reHeaderKeywords = re.compile(r"([-+]?\d+)\s+polygons|pixelization\s+([-+]?\d+)([sd])|balkanized|snapped|real\s+([-+]?\d+)")
        #
        ff = open(filename,"r")
        self.npoly = None
        self.pixelization=None
        self.snapped=False
        self.balkanized=False
        self.real=8
        self.names=[]
        self.metadata={}
        line = ff.readline()
        #loop to read through header lines - stops when rePoly.match(line) returns something, indicating that we've
        #made it through the header and have reached a polygon.
        while (rePoly.search(line) is None)&(len(line)>0):
            sss=reHeaderKeywords.search(line)
            #if line starts with a header keyword, add the info as metadata
            if sss is not None:
                if rePolycount.search(line) is not None:
                    self.npoly = int( sss.group(1) )
                elif reReal.search(line) is not None:
                    self.real = int( sss.group(4) )
                elif rePixelization.search(line) is not None:
                    self.pixelization = (int(sss.group(2)),sss.group(3))
                elif reSnapped.search(line) is not None: 
                    self.snapped = True
                elif reBalkanized.search(line) is not None:
                    self.balkanized = True
            #print warning if the line doesn't match any of the header keywords
            else:
                print "WARNING: line \""+line+"\" in "+filename+" ignored."
            #read next line
            line = ff.readline()
                
        if self.npoly is None:
            raise RuntimeError,"Did not find polygon count line \"n polygons\" in header of %s"%filename
        #if a value for real is given in the keyword args, force it to override what's in the file
        if real is not None:
            self.real=real
        #if self.real = 10 and longdouble_utils imported successfully, use long doubles. Otherwise just use float64
        if self.real==10:
            self.uselongdoubles=uselongdoubles
        elif self.real==8:
            self.uselongdoubles=False
        else:
            raise RuntimeError,"bad value for real: %d.  Should be real=8 to use doubles, real=10 to use long doubles."%self.real
        if self.uselongdoubles:
            floattype=np.longdouble
            self.real=10
        else:
            floattype=np.float64
            self.real=8
        self.polylist = empty(self.npoly,dtype='object')
        self.polyids = zeros(self.npoly,dtype=int)
        self.areas = -ones(self.npoly,dtype=floattype)
        self.weights = zeros(self.npoly,dtype=floattype)
        self.ncaps = zeros(self.npoly,dtype=int)
        self.pixels = zeros(self.npoly,dtype=int)
        counter = 0
        ss = rePoly.search(line)
        while len(line)>0:
            while (ss==None)&(len(line)>0):
                line = ff.readline()
                ss   = rePoly.search(line)
            if len(line)>0:
                ipoly= int(ss.group(1))
                ncap = int(ss.group(2))
                # Check to see if we have a weight.
                ss = reWeight.search(line)
                if ss==None:
                    weight=floattype(0.0)
                else:
                    if self.uselongdoubles:
                        weight=longdouble_utils.string2longdouble(ss.group(1))
                    else:
                        weight=float(ss.group(1))
                # Check to see if we have an area.
                ss = reArea.search(line)
                if ss==None:
                    area= floattype(-1.0)
                else:
                    if self.uselongdoubles:
                        area=longdouble_utils.string2longdouble(ss.group(1))
                    else:
                        area=float(ss.group(1))
                # Check to see if we have a pixel number.
                ss = rePixel.search(line)
                if ss==None:
                    pixel = 0
                else:
                    pixel = float(ss.group(1))
                self.polyids[counter] = ipoly
                self.areas[counter] = area
                self.weights[counter] = weight
                self.ncaps[counter] = ncap
                self.pixels[counter] = pixel
                # NOTE: Looping over a numpy array appears to be slower
                # than a python list, using a numpy array for polylist slows
                # down polyid(),area(),etc. by ~2x.  But it makes
                # the get_XX() functions cleaner.
                self.polylist[counter] = zeros((ncap,4),dtype=floattype)
                for i in range(ncap):
                    line = ff.readline()
                    if self.uselongdoubles:
                        cap = [longdouble_utils.string2longdouble(x) for x in string.split(line)]
                    else:
                        cap  = [float(x) for x in string.split(line)]
                    self.polylist[counter][i] = cap
                ss = None
                counter += 1
        w, = np.where(self.pixels != 0)
        #self.npixels = len(set(self.pixels))
        self.npixels = np.unique(w).size
        ff.close()

        #read in extra columns stored in this format:
        # polygon file named poly.pol or poly_obstime.pol (to denote that the polygons are weighted by observation time)
        # additional columns can be provided as poly.maglims, poly.obstime, poly.XXX where XXX is whatever tag applies to the column
        # an additional column will be added to the mangle polygon object named XXX 
        # additional column files should have number of rows equal to the number of polygons in poly.pol
        # this allows alternate weights and extra information about the polygons that are stored in a fits file to be written in ascii format      
        if read_extra_columns:
            #grab base name of file, e.g., poly.pol will yield base_filename=poly
            base_filename=string.split(self.filename,'.')[0]
            #get list of files that satisfy base_filename.* 
            files=glob.glob(base_filename+'.*')
            #remove any files that are clearly not meant to be extra column files from the list - anything that ends in .pol,.ply,.fits,.list,.list.weight, or .eps
            files=[f for f in files if (f[-4:] != '.ply') & (f[-4:] != '.pol') & (f[-5:] != '.fits') & (f[-5:] != '.list') & (f[-5:] != '.eps') &  (f[-12:] != '.list.weight')]
            #if that didn't find any files to read in, try stripping off everything after the final underscore for the basename, so poly_obstime.pol yields base_filename=poly
            if len(files)==0:
                weights_tag=string.split(base_filename,'_')[-1] 
                base_filename='_'.join(string.split(base_filename,'_')[:-1])
                files=glob.glob(base_filename+'.*')
                files=[f for f in files if (f[-4:] != '.ply') & (f[-4:] != '.pol') & (f[-5:] != '.fits') & (f[-5:] != '.list') & (f[-5:] != '.eps') &  (f[-12:] != '.list.weight')]
            #loop through files
            for f in files:
                try:
                    #basic file read that should work for anything with a fixed number of columns
                    data=np.genfromtxt(f,dtype=None,invalid_raise=True)
                except:
                    try:
                        #this can read a file that has a variable number of integers on each line - e.g., the output of the 'polyid' mangle function
                        data=array([array([int(x) for x in string.split(line)]) for line in open(f,'r')])                        
                    except:
                        try:
                            #this can read a file with a variable number of floats on each line
                            data=array([array([float(x) for x in string.split(line)]) for line in open(f,'r')])
                        except:
                            print 'WARNING: could not read column from file '+f
                            continue                            
                if len(data)!=self.npoly:
                    print 'number of lines in '+f+' does not match number of polygons in '+self.filename
                    continue
                name=string.split(f,'.')[-1]
                self.add_column(name,data)
                
    #...

    def convert_db(self, windows):
        """Read in polygons from a windows db table in photofielddb."""
        self.npoly = len(windows)
        self.pixelization=None

        # pull out the relevant fields
        self.polylist= empty(self.npoly,dtype='object')
        self.polyids= arange(0,self.npoly,dtype=int)
        self.npixels = 0
        self.areas= zeros(self.npoly, dtype= 'float64')
        self.weight= zeros(self.npoly, dtype= 'float32')
        self.pixels= zeros(self.npoly, dtype= 'int32')
        self.ifield= zeros(self.npoly, dtype= 'int32')
        self.ncaps= zeros(self.npoly, dtype='int32')
        for i in range(self.npoly):
            self.areas[i]= windows[i].str
            self.weight[i]= windows[i].weight
            self.pixels[i]= windows[i].pixel
            self.ifield[i]= windows[i].ifield
            self.ncaps[i] = windows[i].ncaps
        for i,n,w in zip(self.polyids,self.ncaps,windows):
            self.polylist[i] = zeros((n,4))
            self.polylist[i][...,:-1] = array(w.xcaps[:3*n]).reshape(n,3)
            self.polylist[i][...,-1] = w.cmcaps[:n]
    #...

    def read_fits_file(self,filename,real=None):
        """Read in polygons from a .fits file."""
        data = pyfits.open(filename,memmap=True)[1].data
        header = pyfits.open(filename,memmap=True)[0].header
        self.npoly = len(data)
        self.pixelization=None
        self.snapped=False
        self.balkanized=False
        self.real=8
        self.names=[]
        self.formats=()
        names = data.dtype.names
        formats=data.formats

        #if mangle header keywords are present, use them
        if ('PIXRES' in header.keys()) & ('PIXTYPE' in header.keys()) &  ('SNAPPED' in header.keys()) &  ('BLKNIZED' in header.keys())  &  ('REAL' in header.keys()):
            self.pixelization=(header['PIXRES'], header['PIXTYPE'])
            self.snapped=header['SNAPPED']
            self.balkanized=header['BLKNIZED']
            self.real=header['REAL']
        #if mangle header keywords aren't present, check for text in header that looks like the header of an ascii polygon file:
        else:
            rePixelization = re.compile(r"pixelization\s+([-+]?\d+)([sd])")
            reSnapped = re.compile(r"snapped")
            reBalkanized = re.compile(r"balkanized")
            reReal = re.compile(r"real\s+([-+]?\d+)")
            reHeaderKeywords = re.compile(r"([-+]?\d+)\s+polygons|pixelization\s+([-+]?\d+)([sd])|balkanized|snapped|real\s+([-+]?\d+)")
            cardlist = header.ascard

            for card in cardlist:
                line=str(card)
                sss=reHeaderKeywords.search(line)
                #if line starts with a header keyword, add the info as metadata
                if sss is not None:
                    if rePixelization.search(line) is not None:
                        result = rePixelization.search(line)
                        self.pixelization = (int(result.group(1)),result.group(2))
                    elif reReal.search(line) is not None:
                        result = reReal.search(line)
                        self.real = int(resulth.group(3))
                    elif reSnapped.search(line) is not None: 
                        self.snapped = True
                    elif reBalkanized.search(line) is not None:
                        self.balkanized = True

        #if a value for real is given in the keyword args, force it to override what's in the file
        if real is not None:
            self.real=real
        #if self.real = 10 and longdouble_utils imported successfully, use long doubles. Otherwise just use float64
        if self.real==10:
            self.uselongdoubles=uselongdoubles
        elif self.real==8:
            self.uselongdoubles=False
        else:
            raise RuntimeError,"bad value for real: %d.  Should be real=8 to use doubles, real=10 to use long doubles."%self.real
        if self.uselongdoubles:
            floattype=np.longdouble
            self.real=10
        else:
            floattype=np.float64
            self.real=8

        # pull out the relevant fields
        self.polylist = empty(self.npoly,dtype='object')
        self.polyids = arange(0,self.npoly,dtype=int)
        if 'STR' in names:
            if ('STR_X' in names) & (self.uselongdoubles):
                self.areas=longdouble_utils.doubledouble2longdouble((data['STR'],data['STR_X']))
            else:
                self.areas = floattype(data['STR'])
        else:
            self.areas = -ones(self.npoly,dtype=floattype)
        if 'WEIGHT' in names:
            if ('WEIGHT_X' in names) & (self.uselongdoubles):
                self.weights=longdouble_utils.doubledouble2longdouble((data['WEIGHT'],data['WEIGHT_X']))
            else:
                self.weights = floattype(data['WEIGHT'])
        else:
            self.weights = zeros(self.npoly,dtype=floattype)
        if 'PIXEL' in names:
            self.pixels = data['PIXEL']
        else:
            self.pixels = zeros(self.npoly,dtype=int)
        self.npixels = len(set(self.pixels))
        self.ncaps = data['NCAPS']
        for i,n,x in izip(self.polyids,self.ncaps,data):
            self.polylist[i] = zeros((n,4),dtype=floattype)
            if ('XCAPS_X' in names) & (self.uselongdoubles):
                xcaps=longdouble_utils.doubledouble2longdouble((x['XCAPS'],x['XCAPS_X']))
            else:
                xcaps=floattype(x['XCAPS'])
            if ('CMCAPS_X' in names) & (self.uselongdoubles):
                cmcaps=longdouble_utils.doubledouble2longdouble((x['CMCAPS'],x['CMCAPS_X']))
            else:
                cmcaps=floattype(x['CMCAPS'])            
            self.polylist[i][...,:-1] = xcaps.reshape(-1,1)[:3*n].reshape(n,3)
            self.polylist[i][...,-1] = cmcaps.reshape(-1)[:n]

        # Read any additional fields that may be in the file
        info=data.columns.info(output=False)
        self.metadata={}
        for i, name in enumerate(names):
            if ((name != 'XCAPS') & (name != 'CMCAPS') & (name != 'NCAPS') & (name != 'STR') & (name != 'WEIGHT') & (name != 'PIXEL')
                & (name != 'XCAPS_X') & (name != 'CMCAPS_X') & (name != 'STR_X') & (name != 'WEIGHT_X')):
               # fits files use bzero and bscale along with signed int data types to represent unsigned integers
                # e.g. for an unsigned 32 bit integer, the format code will be 'J' (same as for a signed 32 bit int),
                # bscale will be 1 and bzero will be 2**31.  The data gets automatically converted as data=bscale*rawdata+bzero
                if (info['bscale'][i] == 1.0) & (info['bzero'][i] is not ''):
                    precision=str(np.int(np.log2(np.abs(np.double(info['bzero'][i])))+1)) # equal to '8', '16','32', '64'
                    if ((precision == '16') | (precision == '32') | (precision == '64')) & (info['bzero'][i]>0):
                        #define a numpy data type that is an unsigned int of the given precision
                        dt=np.typeDict['uint'+precision]
                        col = dt(data[name])
                    #for 8 bit ints, it's vice-versa: fits natively has an unsigned 8 bit int and requires bzero to convert to signed.
                    elif (precision == '8') & (info['bzero'][i]<0):
                        #define a numpy data type that is an unsigned int of the given precision 
                        dt=np.typeDict['int'+precision]
                        col = dt(data[name])
                    else:
                        col = data[name]
                else:
                    col = data[name]
                self.add_column(string.lower(name), col, format=info['format'][i], asciiformat=None, unit=info['unit'][i], null=info['null'][i], bscale=info['bscale'][i], bzero=info['bzero'][i], disp=info['disp'][i], start=info['start'][i], dim=info['dim'][i])


    def add_column(self,name, data, format=None, asciiformat=None, unit=None, null=None, bscale=None, bzero=None, disp=None, start=None, dim=None):
        """
        Add a column to the polygons object. Keyword arguments are the
        same as the pyfits.Column constructor.
        Format will be detected automatically from the data type if not provided.
        """
        #check to make sure length of array matches number of polygons
        if len(data) != self.npoly:
            raise RuntimeError('Length of array does not match number of polygons.')
        #detect dimensions of array
        if data.ndim==1:
            size=''
        elif data.ndim==2:
            size=str(data.shape[-1])
        elif data.ndim>2:
            size=str(np.prod(data.shape[1:]))
            if dim is None:
                dim=str(data.shape[1:])

        if data.dtype.type is np.object_:
            size='P'
            dt=data[0].dtype.type
        else:
            dt=data.dtype.type
                
        #detect appropriate format based on numpy type
        #use bzero and bscale for unsigned integers                  
        if dt is np.float64:
            ftype='D'
            asciiformat='%.15g'
        elif dt is np.int32:
            ftype='J'
            asciiformat='%d'
        elif dt is np.uint32:
            ftype='J'
            bzero=2**31
            bscale=1
            asciiformat='%u'
        elif dt is np.string_:
            #add the string length to the dimensions in order to store an array of strings
            dims=(data.dtype.itemsize,)+data.shape[1:]
            size=str(np.prod(dims))
            if len(dims)>1:
                if dim is None:
                    dim=str(dims)
                ftype='A'
            asciiformat='%s'
        elif dt is np.float32:
            ftype='E'
            asciiformat='%.6g'
        elif dt is np.int16:
            ftype='I'
            asciiformat='%d'
        elif dt is np.uint16:
            ftype='I'
            bzero=2**15
            bscale=1
            asciiformat='%u'
        elif dt is np.int64:
            ftype='K'
            asciiformat='%d'
        elif dt is np.uint64:
            #64 bit unsigned ints don't work right now b/c scaling is done by converting bzero to a double, which doesn't have enough precision
            #format='K'
            #bzero=2**63
            #use 32 bit unsigned int instead:
            type='J'
            bzero=2**31
            bscale=1
            asciiformat='%d'
        elif dt is np.complex64:
            ftype='C'
            asciiformat=None
        elif dt is np.complex128:
            ftype='M'
            asciiformat=None
        elif dt is np.bool_:
            ftype='L'
            asciiformat='%d'
        elif dt is np.int8:
            ftype='B'
            bzero=-2**7
            bscale=1
            asciiformat='%d'
        elif dt is np.uint8:
            ftype='B'
            asciiformat='%u'
        else:
            raise RuntimeError('data type '+str(dt)+' not supported')

        if size=='P':
            asciiformat=None
        if format is None:
            format=size+ftype
            if size=='P':
                format=format+'()'
        
        self.metadata.update({name:{'bscale':bscale,
                                    'bzero':bzero,
                                    'dim':dim,
                                    'disp':disp,
                                    'format':format,
                                    'asciiformat':asciiformat,
                                    'name':string.upper(name),
                                    'null':null,
                                    'start':start,
                                    'unit':unit}})
        vars(self)[name] = data
        self.names+=[name]

    def remove_column(self,name):
        """removes a column from the polygon structure and metadata"""
        del vars(self)[name]
        del self.metadata[name]
        self.names.remove(name)
        
    def write_fits_file(self,filename,clobber=True,keep_ids=False):
        """Write polygons to a .fits file.

        If clobber=True, overwrite existing file with the same name.

        If keep_ids=True, write id numbers in self.ids (or self.polyids
        if self.ids does not exist) as 'ids' column in fits file.        
        """

        #define size of xcaps and cmcaps arrays based on the maximum number of caps in any polygon, fill extra spaces with zeros
        maxn=self.ncaps.max()

        if self.uselongdoubles:
            floattype=np.longdouble
        else:
            floattype=np.float64

        #initialize to zero-filled arrays
        xcaps=zeros((len(self.polylist),maxn,3),dtype=floattype)
        cmcaps=zeros((len(self.polylist),maxn),dtype=floattype)

        #run through polygons and massage caps data in polylist to fit into cmcaps and xcaps arrays
        polyids=arange(0,self.npoly,dtype=int)
        for i,n,poly in izip(polyids,self.ncaps,self.polylist):
            cmcaps[i,:n]= poly[...,-1]
            xcaps[i,:n]=poly[...,:-1]

        #reshape xcaps array to have dimensions (npolys, 3*maxn) rather than (npolys,maxn,3) to keep pyfits happy  
        xcaps=xcaps.reshape(-1,3*maxn)

        if self.uselongdoubles:
            cmcaps_doubdoub=longdouble_utils.longdouble2doubledouble(cmcaps)
            xcaps_doubdoub=longdouble_utils.longdouble2doubledouble(xcaps)
            weights_doubdoub=longdouble_utils.longdouble2doubledouble(self.weights)
            areas_doubdoub=longdouble_utils.longdouble2doubledouble(self.areas)
            cmcaps=cmcaps_doubdoub[0]
            cmcaps_x=cmcaps_doubdoub[1]
            xcaps=xcaps_doubdoub[0]
            xcaps_x=xcaps_doubdoub[1]
            weights=weights_doubdoub[0]
            weights_x=weights_doubdoub[1]
            areas=areas_doubdoub[0]
            areas_x=areas_doubdoub[1]
        else:
            weights=self.weights
            areas=self.areas
            
        #if keeping polygon id number in fits file, add 'ids' column
        if keep_ids == True:
            if 'ids' not in self.names:
                self.add_column('ids',self.polyids)

        #define pyfits columns for basic polygon file elements
        xcaps_col=pyfits.Column(name='XCAPS',format=str(3*maxn)+'D',dim='( 3, '+str(maxn)+')',array=xcaps)
        cmcaps_col=pyfits.Column(name='CMCAPS',format=str(maxn)+'D',array=cmcaps)
        ncaps_col=pyfits.Column(name='NCAPS',format='J',array=self.ncaps)
        str_col=pyfits.Column(name='STR',format='D',array=areas)
        weight_col=pyfits.Column(name='WEIGHT',format='D',array=weights)
        pixel_col=pyfits.Column(name='PIXEL',format='J',array=self.pixels)
        if self.uselongdoubles:
            xcaps_x_col=pyfits.Column(name='XCAPS_X',format=str(3*maxn)+'D',dim='( 3, '+str(maxn)+')',array=xcaps_x)
            cmcaps_x_col=pyfits.Column(name='CMCAPS_X',format=str(maxn)+'D',array=cmcaps_x)
            str_x_col=pyfits.Column(name='STR_X',format='D',array=areas_x)
            weight_x_col=pyfits.Column(name='WEIGHT_X',format='D',array=weights_x)

        #define pyfits columns for any extra columns listed in self.names
        extracols=[]
        for name in self.names:
            info=self.metadata[name]
            data=vars(self)[name]
            #reshape arrays into 1-dim arrays so pyfits can deal
            if data.dtype.type is np.string_:
                if data.ndim>1:
                    raise RuntimeError('arrays of strings not supported for exporting to fits file.')
                   # using the line below give "operands could not be broadcast together" error
                   #  data=array([''.join([a.ljust(data.dtype.itemsize) for a in x.flat]) for x in data])                 
            else:
                if data.ndim>2:
                    data=data.reshape(-1,np.prod(data.shape[1:]))
            col=pyfits.Column(name=info['name'],
                              format=info['format'],
                              bscale=info['bscale'],
                              bzero=info['bzero'],
                              dim=info['dim'],
                              disp=info['disp'],
                              null=info['null'],
                              start=info['start'],
                              unit=info['unit'],
                              array=data)
            extracols=extracols+[col]

        #make a primary HDU for the table and add metadata
        primary=pyfits.PrimaryHDU()
        if self.pixelization is not None:
            primary.header.update('PIXRES',self.pixelization[0])
            primary.header.update('PIXTYPE',self.pixelization[1])    
        primary.header.update('SNAPPED',self.snapped)
        primary.header.update('BLKNIZED',self.balkanized)
        primary.header.update('REAL',self.real)

        #make new fits table and write it to file
        if self.uselongdoubles:
            table=pyfits.new_table(pyfits.ColDefs([xcaps_col,xcaps_x_col,cmcaps_col,cmcaps_x_col,ncaps_col,weight_col,weight_x_col,pixel_col,str_col,str_x_col]+extracols))
        else:
            table=pyfits.new_table(pyfits.ColDefs([xcaps_col,cmcaps_col,ncaps_col,weight_col,pixel_col,str_col]+extracols))
        tablehdus=pyfits.HDUList([primary, table])
        tablehdus.writeto(filename,clobber=clobber)

    def __init__(self,filename,db=False,keep_ids=False,read_extra_columns=False,real=None):
        """
        Initialize Mangle with a file containing the polygon mask.
        If db == True, filename is expected to be a windows db table.

        If keep_ids=True, read in input id numbers from file and keep as
        self.ids column, and use self.polyids as a 0 through npoly index

        If read_extra_columns=True, supplemental files for ascii polygon (.ply/.pol) files
        containing additional information will be read in this format:
        For a polygon file named poly.pol or poly_obstime.pol (to denote that the polygons
        are weighted by observation time), additional columns can be provided as poly.maglims,
        poly.obstime, poly.XXX where XXX is whatever tag applies to the column.
        An additional column will be added to the mangle polygon object named XXX. 
        Additional column files should have number of rows equal to the number of
        polygons in poly.pol. This allows alternate weights and extra information about the
        polygons that are stored in a fits file to be written in ascii format    

        Acceptable formats (determined from the file extension):
            .ply or .pol <-- Mangle polygon files:
                http://space.mit.edu/~molly/mangle/manual/polygon.html
            .fits <-- FITS file
            windows db table <-- if db == True
        """
        if db == True:
            self.convert_db(filename)
            self.filename == None
        else:
            if not os.path.exists(filename):
                raise IOError,"Can not find %s"%filename
            self.filename = filename # useful to keep this around.
            if (filename[-4:] == '.ply') | (filename[-4:] == '.pol'):
                self.read_ply_file(filename,read_extra_columns=read_extra_columns,real=real)
            elif filename[-5:] == '.fits':
                self.read_fits_file(filename,real=real)
            else:
                raise IOError,"Unknown file extension for %s"%filename

        if self.pixelization:
            self.pixel_dict = self._create_pixel_dict()

        if len(self.polylist) != self.npoly:
            print "Got %d polygons, expecting %d."%\
              (len(self.polylist),self.npoly)
        if keep_ids == True:
            if 'ids' not in self.names:
                self.add_column('ids',self.polyids)
                self.polyids=arange(0,self.npoly,dtype=int)
        else:
            # Check whether the polyids are sequential and range from 0 to npoly-1
            # If they don't, then there may be a problem with the file.
            # NOTE: this should always be correct for current_boss_geometry.
            badcounter = 0
            if (min(self.polyids)==0) & (max(self.polyids)==self.npoly-1) :
                for i,polyid in enumerate(self.polyids):
                    if i != polyid:
                        badcounter += 1
                if badcounter > 0:
                    print "WARNING!!!!"
                    print "Found",badcounter,"polygons out of order."
                    print "Reordering polygons so that polyid=index"
                    self.sortpolys()
                    badcounter = 0
            else:
                print "WARNING!!!!"
                print "Range of polygon ids in input is (",min(self.polyids),",",max(self.polyids),"), not ( 0 ,",self.npoly-1,")"
                print "Forcing 'polyids' attribute to be 0 to npoly-1 and saving ids from input as 'id' attribute"
                print "To do this automatically, use 'keep_ids=True' in the mangle.Mangle constructor."
                print "To write polygon file retaining the input ids, use 'keep_ids=True' in writeply()."
                self.add_column('ids',self.polyids)
                self.polyids=arange(0,self.npoly,dtype=int)
    #...
     
    def _create_pixel_dict(self):
        """
        Return dictionary of what indexes are in each pixel,
        to improve pixel-based polyid searches.
        """
        start = time.time()
        pixel_dict = {}
        for i,pix in enumerate(self.pixels):
            if pix not in pixel_dict:
                pixel_dict[pix] = []
            pixel_dict[pix].append(i)
        end = time.time()
        #print "!!!pixel_dict creation time:",end-start
        return pixel_dict
     #...
#...
