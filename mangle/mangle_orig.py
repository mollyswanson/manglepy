#!/usr/bin/env python
# Python class to implement some basic Mangle routines
# to handle caps etc.
# Stand-alone python, no non-standard Python modules required.
#
# This routine is slow, to work with large amounts of data
# use compiled routines instead.
#
# Author:		Martin White	(UCB)
# Written		17-Mar-2010
# Modified:		17-Mar-2010	(Basic bug fixes)
#			22-Mar-2010	(Track pixel area)
#			20-Apr-2010	(Add setweight and writeply methods)
#			21-Apr-2010	(Add totalarea method)
#			12-May-2010	(Allow no space before "caps")
#

import os
import re
import string
import math as M


class Mangle:

    """
    Mangle:
    A Python class to implement some basic mangle routines to check
    what polygon a given (ra,dec) lies within, what its Mangle weight
    is etc.  [(ra,dec) should be in decimal degrees.]
    The class is initialized with 1 argument, the file name of an
    ascii string (in Mangle polygon format) containing the mask.
    If the file contains no weights, they are set to zero for all polygons.
    """

    __author__ = "Martin White"
    __version__ = "1.0"
    __email__  = "mwhite@berkeley.edu"

    def incap_spam(self,cap,theta,phi):
        """
        incap_spam(self,cap,theta,phi):
        This is an internal function used by the class, you shouldn't need
        to use it.
        Returns True if (theta,phi) lies in the cap specified by the
        4-vector "cap" containing (x,y,z,cm).
        """
        x0 = M.sin(theta)*M.cos(phi)
        y0 = M.sin(theta)*M.sin(phi)
        z0 = M.cos(theta)
        cd = 1.0-cap[0]*x0-cap[1]*y0-cap[2]*z0
        if cap[3]<0.0:
            if cd>M.fabs(cap[3]):
                tmp = True
            else:
                tmp = False
        else:
            if cd<cap[3]:
                tmp = True
            else:
                tmp = False
        return(tmp)

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
        tmp = True
        for cap in polygon[3:]:
            tmp = tmp&self.incap_spam(cap,theta,phi)
        return(tmp)

    def polyid(self,ra,dec):
        """
        polyid(self,ra,dec):
        This is one of the main methods for this class, it returns the
        ID number of the polygon given an RA and DEC in decimal degrees.
        """
        theta = M.pi/180. * (90.0-dec)
        phi   = M.pi/180. * ra
        ipoly = -1
        for poly in self.polylist:
            if self.inpoly_spam(poly,theta,phi):
                ipoly = poly[0]
                break
        return(ipoly)

    def weight(self,ra,dec):
        """
        weight(self,ra,dec):
        This is one of the main methods for this class, it returns the
        polygon weight given an RA and DEC in decimal degrees.
        """
        theta = M.pi/180. * (90.0-dec)
        phi   = M.pi/180. * ra
        weight= -1
        for poly in self.polylist:
            if self.inpoly_spam(poly,theta,phi):
                weight = poly[1]
                break
        return(weight)

    def area(self,ra,dec):
        """
        area(self,ra,dec):
        This is one of the main methods for this class, it returns the
        polygon area given an RA and DEC in decimal degrees.
        """
        theta = M.pi/180. * (90.0-dec)
        phi   = M.pi/180. * ra
        area  = -1.0
        for poly in self.polylist:
            if self.inpoly_spam(poly,theta,phi):
                area = poly[2]
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
        for poly in self.polylist:
            tot_area += poly[2]
            eff_area += poly[2]*poly[1]
        return((tot_area,eff_area))

    def setweight(self,polyid,weight):
        """
        setweight(self,polyid,weight):
        Sets the weight entry of polygon "polyid" to weight.
        """
        for i,poly in enumerate(self.polylist):
            if poly[0]==polyid:
                self.polylist[i][1]=weight
                break

    def setallweight(self,weight):
        """
        setallweight(self,weight):
        Sets the weight entry of all polygons to weight, typically used
        to set everything to 0.
        """
        for i in range(len(self.polylist)):
            self.polylist[i][1]=weight


    def writeply(self,fn):
        """
        writeply(self,fn):
        Writes a Mangle-formatted polygon file containing the information
        in the class.
        """
        ff = open(fn,"w")
        ff.write("%d polygons\n"%len(self.polylist))
        for poly in self.polylist:
            str = "polygon %10d ( %d caps,"%(poly[0],len(poly)-3)
            str+= " %.8f weight, %.15f str):\n"%(poly[1],poly[2])
            ff.write(str)
            for cap in poly[3:]:
                ff.write("%25.20f %25.20f %25.20f %25.20f\n"%\
                  (cap[0],cap[1],cap[2],cap[3]))
        ff.close()

    def __init__(self,fn):
        """
        __init__(self,fn):
        The class is initialized with the name of an ascii file containing
        the Mangle mask.
        """
        if not os.path.exists(fn):
            raise RuntimeError,"Can not find %s"%fn
        #
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
        #
        self.polylist = []
        #
        ss = ex1.match(line)
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
                polyg= [ipoly,weight,area]
                for i in range(ncap):
                    line = ff.readline()
                    cap  = [float(x) for x in string.split(line)]
                    polyg.append(cap)
                self.polylist.append(polyg)
                ss=None
        ff.close()
        if len(self.polylist) != self.npoly:
            print "Got %d polygons, expecting %d."%\
              (len(self.polylist),self.npoly)
        #
