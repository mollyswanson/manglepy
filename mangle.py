#!/usr/bin/env python
"""
Python class to implement some basic Mangle routines
to handle caps etc.

The single element versions of this are fairly slow when looping over
a large data set, but the get_polyids/get_areas functions are within
a factor of two of the speed of the commandline code.

Requires numpy, and pyfits > 2.3.1 for loading fits files.
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

import os
import re
import string

import numpy as np
# some of these will need to be called often, so explicitly name them
from numpy import array,sin,cos,fabs,pi,empty,zeros,ones,argsort,take,arange

class Mangle:
    """
    Mangle:
    A Python class to implement some basic mangle routines to check
    what polygon a given (ra,dec) lies within, what its Mangle weight
    is etc.  [(ra,dec) should be in decimal degrees.]
    The class is initialized with 1 argument, the file name of an
    ascii string (in Mangle polygon format) containing the mask.
    If the file contains no weights, they are set to zero for all polygons.
    
    Warning: the polyid(), area(), and weight() functions are very
    slow when used in a loop.  get_polyids(), etc. are up to 40x
    faster and should be used anytime more than ~100 objects need to
    have their polygons found.
    """

    __author__ = "John Parejko & Martin White"
    __version__ = "2.0"
    __email__  = "john.parejko@yale.edu"

    def incap_spam(self,cap,x0,y0,z0):
        """
        incap_spam(self,cap,x0,y0,z0):
        This is an internal function used by the class, you shouldn't need
        to use it.
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
        Returns True for each (theta,phi) that lies in the cap specified by the
        4-vector "cap" containing (x,y,z,cm), and False for the rest.
        """
        cd = 1.0-cap[0]*x0-cap[1]*y0-cap[2]*z0
        return ((cap[3] < 0.0) & (cd>fabs(cap[3]))) | ((cap[3] > 0.0) & (cd<cap[3]))
    #...

    def inpoly_spam(self,polygon,theta,phi):
        """
        inpoly_spam(self,polygon,theta,phi):
        This is an internal function used by the class, you shouldn't need
        to use it.
        Returns True if (theta,phi) is in the polygon, i.e. if it is
        within all of the caps in the polygon.
        A polygon is a list giving the polygon number, the weight and
        then all of the caps (each cap is a 4-vector/4-list).
        """
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
        inpoly_spam(self,polygon,theta,phi):
            Returns True if (theta,phi) is in the polygon, i.e. if it is
            within all of the caps in the polygon.
            A polygon is a list giving the polygon number, the weight and
            then all of the caps (each cap is a 4-vector/4-list).
        """
        test = ones(len(x0),dtype=bool)
        for x in polygon:
            test &= self.incap_vec(x,x0,y0,z0)
        return test
        # NOTE: the below doesn't work, but I think it should.
        #return all((self.incap_vec(x,x0,y0,z0) for x in polygon[3]),axis=1)
    #...

    def get_polyids(self,ra,dec):
        """
        polyid(self,ra,dec):
            Return the ID number of polygons given arrays of RA and DEC
            in decimal degrees.
        """
        theta = pi/180. * (90.0-dec)
        phi = pi/180. * ra
        sintheta = sin(theta)
        x0 = sintheta*cos(phi)
        y0 = sintheta*sin(phi)
        z0 = cos(theta)
        goodpolys = -ones(len(ra),dtype=int)
        for i,poly in zip(self.polyids,self.polylist):
            test = self.inpoly_vec(poly,x0,y0,z0)
            goodpolys[test] = i
        return goodpolys
    #...

    def get_areas(self,ra,dec):
        """
        get_areas(self,ra,dec):
            Return the areas of the polygons containing each RA/Dec pair.
            Result is in steradians.
        """
        polyids=self.get_polyids(ra,dec)
        return self.areas[polyids]
    #...

    def get_weights(self,ra,dec):
        """
        get_weights(self,ra,dec):
            Return the weights of the polygons containing each RA/dec pair.
        """
        polyids=self.get_polyids(ra,dec)
        return self.weights[polyids]
    #...

    def polyid(self,ra,dec):
        """
        polyid(self,ra,dec):
        This is one of the main methods for this class, it returns the
        ID number of the polygon given an RA and DEC in decimal degrees.
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
        """
        weight(self,ra,dec):
        This is one of the main methods for this class, it returns the
        polygon weight given an RA and DEC in decimal degrees.
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
        """
        area(self,ra,dec):
        This is one of the main methods for this class, it returns the
        polygon area in steradians given an RA and DEC in decimal degrees.
        """
        theta = pi/180. * (90.0-dec)
        phi   = pi/180. * ra
        area  = -1.0
        for a,poly in zip(self.areas,self.polylist):
            if self.inpoly_spam(poly,theta,phi):
                area = a
                break
        return(area)

    def totalarea(self):
        """
        totalarea(self):
        Returns the total area in the mask (i.e. the sum of the areas of
        each polygon) and the total "effective" area (i.e. the area weighted
        by the completeness).
        Returns (tot_area,eff_area).
        """
        tot_area,eff_area = 0.0,0.0
        for a,w in zip(self.areas,self.weights):
            """poly in self.polylist:"""
            tot_area += a
            eff_area += a*w
        return((tot_area,eff_area))

    def setweight(self,polyid,weight):
        """
        setweight(self,polyid,weight):
        Sets the weight entry of polygon "polyid" to weight.
        """
        self.weights[polyid] = weight

    def setallweight(self,weight):
        """
        setallweight(self,weight):
        Sets the weight entry of all polygons to weight, typically used
        to set everything to 0.
        """
        self.weights.flat = weight

    def writeply(self,fn):
        """
        writeply(self,fn):
        Writes a Mangle-formatted polygon file containing the information
        in the class.
        """
        ff = open(fn,"w")
        ff.write("%d polygons\n"%len(self.polylist))
        for i in self.polyids:
            str = "polygon %10d ( %d caps,"%(i,len(self.polylist[i]))
            str+= " %.8f weight, %.15f str):\n"%(self.areas[i],self.weights[i])
            ff.write(str)
            for cap in self.polylist[i]:
                ff.write("%25.20f %25.20f %25.20f %25.20f\n"%\
                  (cap[0],cap[1],cap[2],cap[3]))
        ff.close()

    def __getitem__(self,idx):
        """Return an ndarray containing the nonzero attributes[idx]."""
        import numpy.core.records as rec
        names = ['polyid','area','weight','ncap','cap']
        types = [self.polyids.dtype.str,self.areas.dtype.str,self.weights.dtype.str,self.ncaps.dtype.str,self.polylist.dtype.str]
        temp = [self.polyids[idx],self.areas[idx],self.weights[idx],self.ncaps[idx],self.polylist[idx]]
        if 'sector' in dir(self) and self.sector != []:
            names.append('sector')
            types.append(self.sector.dtype.str)
            temp.append(self.sector[idx])
        if 'mmax' in dir(self) and self.mmax != []:
            names.append('mmax')
            types.append(self.mmax.dtype.str)
            temp.append(self.mmaxr[idx])
        if 'fgotmain' in dir(self) and self.fgotmain != []:
            names.append('fgotmain')
            types.append(self.fgotmain.dtype.str)
            temp.append(self.fgotmain[idx])
        # Have to handle things differently if we are accessing only one item
        # because rec.array doesn't work if the list elements are length 1.
        if type(idx) == int:
            result = rec.recarray(shape=2,formats=types,names=names)
            for i in range(len(temp)):
                # NOTE: BROKEN! This doesn't work correctly for some reason
                # for items beyond temp[0].  I don't know why, but 
                # the recarray gets assigned nonsensical values.
                result[0][i] = temp[i]
            result = result[0]
        else:
            result = rec.array(temp,formats=types,names=names)
        return result
    #...

    def graphics(self):
        """Return an array of edge points that can be used to plot
        these polygons.
        Calls 'poly2poly' to generate a temporary graphics file and
        return the output.
        The returned array has len() == npoly, with each element
        being an array of ra/dec pairs. Plot a polygon with, e.g.:
            polys = mng.graphics()
            pylab.plot(polys[0]['ra'],polys[0]['dec'])
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

    def read_ply_file(self,fn):
        """Read in polygons from a .ply file."""
        # It's useful to pre-compile a regular expression for a mangle line
        # defining a polygon.
        ex1 = re.compile(r"polygon\s+(\d+)\s+\(\s*(\d+)\s+caps")
        ex2 = re.compile(r"(\d*\.?\d+)\s+weight")
        ex3 = re.compile(r"(\d*\.?\d+)\s+str")
        #
        ff = open(fn,"r")
        self.npoly = 0
        line = ff.readline()
        ss = re.match(r"(\d+)\s+polygons",line)
        if ss==None:
            raise RuntimeError,"Can not parse 1st line of %s"%fn
        else:
            self.npoly = int( ss.group(1) )
        self.polylist = empty(self.npoly,dtype='object')

        ss = ex1.match(line)
        self.polyids = zeros(self.npoly,dtype=int)
        self.areas = -ones(self.npoly)
        self.weights = zeros(self.npoly)
        self.ncaps = zeros(self.npoly)
        counter = 0
        while len(line)>0:
            while (ss==None)&(len(line)>0):
                line = ff.readline()
                ss   = ex1.match(line)
            if len(line)>0:
                ipoly= int(ss.group(1))
                ncap = int(ss.group(2))
                # Check to see if we have a weight.
                ss = ex2.search(line)
                if ss==None:
                    weight=0.0
                else:
                    weight=float(ss.group(1))
                # Check to see if we have an area.
                ss = ex3.search(line)
                if ss==None:
                    area= -1.0
                else:
                    area=float(ss.group(1))
                self.polyids[counter] = ipoly
                self.areas[counter] = area
                self.weights[counter] = weight
                self.ncaps[counter] = ncap
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
        ff.close()
    #...

    def read_fits_file(self,fn):
        """Read in polygons from a .fits file."""
        import pyfits
        data = pyfits.open(fn,memmap=True)[1].data
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

        self.ncaps = data['NCAPS']
        for i,n,x in zip(self.polyids,self.ncaps,data):
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

    def __init__(self,fn):
        """
        __init__(self,fn):
        The class is initialized with the name of an ascii file containing
        the Mangle mask.  Acceptable formats are .ply and .fits.
        """
        if not os.path.exists(fn):
            raise RuntimeError,"Can not find %s"%fn
        self.filename = fn # useful to keep this around.
        
        if fn[-4:] == '.ply':
            self.read_ply_file(fn)
        elif fn[-5:] == '.fits':
            self.read_fits_file(fn)
        else:
            raise RuntimeError,"Unknown file extension for %s"%fn

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
