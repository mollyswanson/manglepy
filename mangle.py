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
    polyids = mng.get_areas(ra,dec)
    mng.areas[polyids] == mng.get_areas(ra,dec)
    mng.weights[polyids] == mng.get_weights(ra,dec)
    first_mng = mng[:10] # new mangle instance with just the first 10 polygons

Requires numpy > 1.0, and pyfits > 2.3.1 for loading fits files.
"""
# Author:		Martin White	(UCB)
# Written		17-Mar-2010
# Modified:		17-Mar-2010	(Basic bug fixes)
#			22-Mar-2010	(Track pixel area)
#			20-Apr-2010	(Add setweight and writeply methods)
#			21-Apr-2010	(Add totalarea method)
#			12-May-2010	(Allow no space before "caps")
#        24-Jun-2010 jkp: (tweaks to gain ~40% speed improvement.)
#        29-Jun-2010 jkp: (numpy vectorized polyid for >10x speed improvement.)
#        22-Jul-2010 jkp: (completed and tested numpy version. Speed is ~1/2 that of the commandline mangle)
#        10-Aug-2010 jkp: (Added resorting of out-of-order polygons and altered original polylist structure)
#        10-Oct-2010 jkp: (Added __getitem__, returning a new Mangle instance given an index)
#        10-Mar-2011 jkp: (Added support for mangle_utils, with a 4x faster cythonized _incap)
#        18-May-2011 apw: (Added support using the photofield database tables in dr8db/spectradb)

import os
import re
import string
import copy
from itertools import izip

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
from numpy import array,sin,cos,fabs,pi,empty,zeros,ones,argsort,take,arange
import pyfits

class Mangle:
    """Implements basic mangle routines to identify polygons containing points.
    
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
    when operating on numpy arrays.
    
    Supports slicing, e.g.:
        newmng = mng[:10]
        newmng2 = mng[mng.areas < 1e-7]
        newmng3 = mng[[1,3,5,7,9]]
    """

    __author__ = "John Parejko & Martin White"
    __version__ = "2.2"
    __email__  = "john.parejko@yale.edu"

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
            for cap in polygon:
                test &= mangle_utils._incap(cap,x0,y0,z0)
        else:
            for cap in polygon:
                test &= self.incap_vec(cap,x0,y0,z0)
        return test
    #...
    
    def which_pixel(self,ra,dec):
        """Return the pixel numbers for each pair of ra,dec.
        
        UNFINISHED!!!
        The pixelization information is , given pixelization
        resolution res and scheme 's' or 'd'.
        !!! NOTE: only scheme 's' is currently implemented"""
        if self.pixelization == None:
            raise TypeError('No pixelization defined in this mangle instance.')
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ### TBD!!!
        ### Math for this comes from which_pixel in the mangle2.2/src directory
        return None
    #...        

    def get_polyids(self,ra,dec):
        """Return the ID numbers of the polygons containing each RA/Dec pair.

        ra/dec should be numpy arrays in decimal degrees.
        """
        # force an upconversion to double, in case ra/dec are float32
        theta = pi/180. * (90.0-np.float64(dec))
        phi = pi/180. * np.float64(ra)
        sintheta = sin(theta)
        x0 = sintheta*cos(phi)
        y0 = sintheta*sin(phi)
        z0 = cos(theta)
        # If we have a pixelized mask, we can reduce the number of polygons 
        # we have to check for each object.
        #if self.npixels > 0:
        #    for pix in self.pixels:
        #        pixels = self.which_pixel(ra,dec)
        #    # TBD: This needs to be finished!
        #    return None
        #else:
        goodpolys = -ones(len(ra),dtype=int)
        for i,poly in izip(self.polyids,self.polylist):
            test = self.inpoly_vec(poly,x0,y0,z0)
            goodpolys[test] = i
        return goodpolys
    #...

    def get_areas(self,ra,dec):
        """Return the areas of the polygons containing each RA/Dec pair.

        ra/dec should be numpy arrays in decimal degrees.

        Return value is in steradians.

        Polygon areas are taken directly from the input file.
        """
        polyids=self.get_polyids(ra,dec)
        return self.areas[polyids]
    #...

    def get_weights(self,ra,dec):
        """Return the weights of the polygons containing each RA/dec pair.

        ra,dec should be numpy arrays in decimal degrees.

        The weight is taken directly from the input file."""
        polyids=self.get_polyids(ra,dec)
        return self.weights[polyids]
    #...

    def polyid(self,ra,dec):
        """Return the polyid of the polygon containing (ra,dec).
        
        ra,dec should be in decimal degrees.
        """
        theta = pi/180. * (90.0-dec)
        phi   = pi/180. * ra
        ipoly = -1
        for i,poly in zip(self.polyids,self.polylist):
            if self.inpoly_spam(poly,theta,phi):
                ipoly = i
                break
        return(ipoly)

    def weight(self,ra,dec):
        """Return the weight of the polygon containing (ra,dec).
        
        ra,dec should be in decimal degrees.
        """
        theta = pi/180. * (90.0-dec)
        phi   = pi/180. * ra
        weight= -1
        for w,poly in zip(self.weights,self.polylist):
            if self.inpoly_spam(poly,theta,phi):
                weight = w
                break
        return(weight)

    def area(self,ra,dec):
        """Return the area of the polygon containing (ra,dec).
        
        ra,dec should be in decimal degrees.

        Return value is in steradians.
        """
        theta = pi/180. * (90.0-dec)
        phi   = pi/180. * ra
        area  = -1.0
        for a,poly in izip(self.areas,self.polylist):
            if self.inpoly_spam(poly,theta,phi):
                area = a
                break
        return(area)

    def totalarea(self):
        """Return the total area and the 'effective' area in the mask.

        total area is the sum of the areas of each polygon.
        effective area is the area weighted by the completeness.
        """
        tot_area = self.areas.sum()
        eff_area = (self.weights*self.areas).sum()
        return((tot_area,eff_area))

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
        """Set the area of all polygons to weight."""
        self.areas.flat = area

    def set_allweight(self,weight):
        """Set the weight of all polygons to weight."""
        self.weights.flat = weight

    def set_allpixel(self,pixel):
        """Set the weight of all polygons to weight."""
        self.pixels.flat = pixel

    def writeply(self,filename):
        """Write a Mangle-formatted polygon file containing these polygons."""
        ff = open(filename,"w")
        ff.write("%d polygons\n"%len(self.polylist))
        for i in range(self.npoly):
            # TBD: add writing pixel information to the output file.
            str = "polygon %10d ( %d caps,"%(self.polyids[i],len(self.polylist[i]))
            str+= " %.8f weight, %.15f str):\n"%(self.weights[i],self.areas[i])
            ff.write(str)
            for cap in self.polylist[i]:
                ff.write("%25.20f %25.20f %25.20f %25.20f\n"%\
                  (cap[0],cap[1],cap[2],cap[3]))
        ff.close()

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
        if 'sector' in dir(mng2) and mng2.sector != []:
            mng2.sector = mng2.sector[idx]
        if 'mmax' in dir(mng2) and mng2.mmax != []:
            mng2.mmax = mng2.mmax[idx]
        if 'fgotmain' in dir(mng2) and mng2.fgotmain != []:
            mng2.fgotmain = mng2.fgotmain[idx]
        return mng2
    #...

    def graphics(self):
        """Return an array of edge points, to plot these polygons.

        Calls the command-line 'poly2poly' to generate a temporary graphics
        file and return the output.

        The returned array has len() == npoly, with each element
        being an array of ra/dec pairs. Plot a polygon with, e.g.:
            polys = mng.graphics()
            pylab.plot(polys[0]['ra'],polys[0]['dec'])
        or all of them (this may take a while,if there are a lot of polygons):
            polys = mng.graphics()
            for poly in polys:
                pylab.plot(poly['ra'],poly['dec'])
        """
        import tempfile
        import subprocess
        tempIn = tempfile.NamedTemporaryFile('rw+b')
        self.writeply(tempIn.name)
        tempOut = tempfile.NamedTemporaryFile('rw+b')
        call = ' '.join(('poly2poly -q -ol30',tempIn.name,tempOut.name))
        subprocess.call(call,shell=True)
        # NOTE: poly2poly -ol always outputs a separate weight file.
        os.remove(tempOut.name+'.weight')
        tempIn.close()
        # the mangle list format has polygons delimited by NaN pairs.
        dtype=[('ra','>f8'),('dec','>f8')]
        data=np.loadtxt(tempOut,dtype=dtype)
        tempOut.close()
        i = 0
        polys = np.empty(self.npoly,dtype='object')
        temp = []
        for x in data:
            if np.isnan(x[0]) and np.isnan(x[1]):
                polys[i] = np.array(temp,dtype=dtype)
                temp = []
                i += 1
            else:
                temp.append(x)
        return polys
    #...

    def read_ply_file(self,filename):
        """Read in polygons from a .ply file."""
        # It's useful to pre-compile a regular expression for a mangle line
        # defining a polygon.
        rePoly = re.compile(r"polygon\s+(\d+)\s+\(\s*(\d+)\s+caps")
        rePixelization = re.compile(r"pixelization\s+(\d+)([sd])")
        reWeight = re.compile(r"(\d*\.?\d+)\s+weight")
        reArea = re.compile(r"(\d*\.?\d+)\s+str")
        rePixel = re.compile(r"(\d*)\s+pixel")
        #
        ff = open(filename,"r")
        self.npoly = 0
        line = ff.readline()
        ss = re.match(r"(\d+)\s+polygons",line)
        if ss==None:
            raise RuntimeError,"Can not parse 1st line of %s"%filename
        else:
            self.npoly = int( ss.group(1) )
        self.polylist = empty(self.npoly,dtype='object')

        self.polyids = zeros(self.npoly,dtype=int)
        self.areas = -ones(self.npoly)
        self.weights = zeros(self.npoly)
        self.ncaps = zeros(self.npoly)
        self.pixels = zeros(self.npoly,dtype=int)
        counter = 0
        # Pixelization information should be on the next line
        line = ff.readline()
        ss = rePixelization.match(line)
        if ss != None:
            self.pixelization = (int(ss.group(1)),ss.group(2))
        else:
            self.pixelization = None
        ss = rePoly.match(line)
        while len(line)>0:
            while (ss==None)&(len(line)>0):
                line = ff.readline()
                ss   = rePoly.match(line)
            if len(line)>0:
                ipoly= int(ss.group(1))
                ncap = int(ss.group(2))
                # Check to see if we have a weight.
                ss = reWeight.search(line)
                if ss==None:
                    weight=0.0
                else:
                    weight=float(ss.group(1))
                # Check to see if we have an area.
                ss = reArea.search(line)
                if ss==None:
                    area= -1.0
                else:
                    area=float(ss.group(1))
                # Check to see if we have a pixel number.
                ss = rePixel.search(line)
                if ss==None:
                    pixel = 0
                else:
                    pixel=float(ss.group(1))
                self.polyids[counter] = ipoly
                self.areas[counter] = area
                self.weights[counter] = weight
                self.ncaps[counter] = ncap
                self.pixels[counter] = pixel
                # NOTE: Looping over a numpy array appears to be slower
                # than a python list, using a numpy array for polylist slows
                # down polyid(),area(),etc. by ~2x.  But it makes
                # the get_XX() functions cleaner.
                self.polylist[counter] = zeros((ncap,4))
                for i in range(ncap):
                    line = ff.readline()
                    cap  = [float(x) for x in string.split(line)]
                    self.polylist[counter][i] = cap
                ss=None
                counter += 1
        self.npixels = len(set(self.pixels))
        ff.close()
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

    def read_fits_file(self,filename):
        """Read in polygons from a .fits file."""
        data = pyfits.open(filename,memmap=True)[1].data
        self.npoly = len(data)

        names = data.dtype.names
        # pull out the relevant fields
        self.polylist = empty(self.npoly,dtype='object')
        self.polyids = arange(0,self.npoly,dtype=int)
        if 'STR' in names:
            self.areas = data['STR']
        else:
            self.areas = -ones(self.npoly)
        if 'WEIGHT' in names:
            self.weights = data['WEIGHT']
        else:
            self.weights = zeros(self.npoly)
        if 'PIXEL' in names:
            self.pixels = data['PIXEL']
        else:
            self.pixels = zeros(self.npoly,dtype=int)
        self.npixels = len(set(self.pixels))

        self.ncaps = data['NCAPS']
        for i,n,x in izip(self.polyids,self.ncaps,data):
            self.polylist[i] = zeros((n,4))
            self.polylist[i][...,:-1] = x['XCAPS'][:3*n].reshape(n,3)
            self.polylist[i][...,-1] = x['CMCAPS'][:n]

        # Some additional fields that may be in the file
        if 'SECTOR' in names:
            self.sector = data['SECTOR']
        else:
            self.sector = []
        if 'MMAX' in names:
            self.mmax = data['MMAX']
        else:
            self.mmax = []
        if 'FGOTMAIN' in names:
            self.fgotmain = data['FGOTMAIN']
        else:
            self.fgotmain = []
    #...

    def __init__(self,filename,db=False):
        """
        Initialize Mangle with a file containing the polygon mask.
        If db == True, filename is expected to be a windows db table.

        Acceptable formats (determined from the file extension):
            .ply <-- Mangle polygon files:
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
            if filename[-4:] == '.ply':
                self.read_ply_file(filename)
            elif filename[-5:] == '.fits':
                self.read_fits_file(filename)
            else:
                raise IOError,"Unknown file extension for %s"%filename

        if len(self.polylist) != self.npoly:
            print "Got %d polygons, expecting %d."%\
              (len(self.polylist),self.npoly)
        # Check whether the polyids are sequential and range from 0 to npoly
        # If they don't, then there may be a problem with the file.
        # NOTE: this should always be correct for current_boss_geometry.
        badcounter = 0
        for i,polyid in enumerate(self.polyids):
            if i != polyid:
                badcounter += 1
        if badcounter > 0:
            print "WARNING!!!!"
            print "Found",badcounter,"polygons out of order."
            print "Reordering polygons so that polyid=index"
            sorter = argsort(self.polyids)
            self.polylist = take(self.polylist,sorter)
            self.weights = self.weights[sorter]
            self.areas = self.areas[sorter]
            self.polyids = self.polyids[sorter]
    #...
#...
