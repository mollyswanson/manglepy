#!/usr/bin/env python
"""
Test class for the mangle version of python.

The goal is to test both the speed and correctness of improvements to
mangle.py, compared to both Martin White's original version, and to the
known-good, very fast C/Fortran version.

Initial Version: 2010.07.14 John K. Parejko
"""

import unittest
import numpy
import os
import pyfits
import csv
import subprocess
import time
import tempfile
import glob

import mangle
import mangle_orig

class TestMangle(unittest.TestCase):
    def setUp(self):
        # set up the mangle polygon geometries
        #self.polyfile = '../data/current_boss_geometry.ply'
        #self.polyfile = '../data/geometry/geometry-boss_1-11-balkanized.ply'
        #self.polyfile = '../data/geometry/geometry-boss_7-10.ply'
        #self.polyfile = '/Users/parejkoj/astronomy/BOSS/data/geometry/geometry-boss_7-11-balkanized.ply'
        self.polyfile = '/Users/parejkoj/astronomy/BOSS/data/geometry/geometry-boss7.ply'
        self.mng = mangle.Mangle(self.polyfile)
        self.mng_orig = mangle_orig.Mangle(self.polyfile)

        # just require that the fits file exist for now.
        self.fitsfile = self.polyfile.split('.ply')[0]+'.fits'
        if glob.glob(self.fitsfile):
            self.mng_fits = mangle.Mangle(self.fitsfile)
        else:
            self.mng_fits = None

        # read in some test data of galaxies
        # these two lines give a relatively quick test
        data = pyfits.open('/Users/parejkoj/astronomy/BOSS/data/final-boss7.fits')
        self.data = data[1].data[:1000]
        # these two lines give a very long test, and are not recommended
        # for use with any non-vectorized mangle, as they will take for ever.
        #data = pyfits.open('/Users/parejkoj/astronomy/BOSS/data/samples/boss+legacy_7-12_CMASS-randoms-2000000.fits')
        #self.data = data[1].data
        data = None

        self.writetestfile = './mangle_test_polys.ply'
        
        # Write out a file containing RA DEC for the commandline mangle
        self.mangletestfile = './mangle_test_data.dat'
        outfile = file(self.mangletestfile,'w')
        outcsv = csv.writer(outfile,delimiter=' ')
        outfile.write('#RA Dec\n')
        outfile.write('#This file written by test_mangle for unittesting.\n')
        outfile.write('#It may be safely deleted.\n')
        outcsv.writerows(zip(self.data.RA,self.data.DEC))
    #...

    def tearDown(self):
        os.remove(self.mangletestfile)
        if glob.glob(self.writetestfile):
            os.remove(self.writetestfile)
    #...

    def do_polyid_cmd(self):
        """Run the commandline version of polyid.
        Return the polygons and the total time the command took."""
        polyidoutfile = './mangle_test_polyid.dat'
        polyid_call = 'polyid -q '+' '.join((self.polyfile,self.mangletestfile,polyidoutfile))
        # NOTE: have to use time.time() here because of clock() returns the
        # CPU time of python, not of external programs.
        start = time.time()
        subprocess.call(polyid_call,shell=True)
        elapsed = (time.time() - start)
        incsv = csv.reader(file(polyidoutfile),delimiter=' ',skipinitialspace=True)
        temp = incsv.next() # strip off the one line header
        polyids = []
        for x in incsv:
            if len(x) == 3:
                polyids.append(float(x[2])) # lines are: RA,Dec,id
            else:
                polyids.append(-1)
        os.remove(polyidoutfile)
        return numpy.array(polyids),elapsed
    #...
    
    def test_totalarea(self):
        start = time.clock()
        total1 = self.mng_orig.totalarea()
        elapsed1 = (time.clock() - start)
        
        start = time.clock()
        total2 = self.mng.totalarea()
        elapsed2 = (time.clock() - start)
        
        self.assertEqual(total1,total2)
        print 'Elapsed time for total area (orig,new):',elapsed1,elapsed2
    #...

    def test_fits_polyid(self):
        """Compare reading data from a .ply file with a .fits file: polyids"""
        ply = self.mng.get_polyids(self.data.RA,self.data.DEC)
        fits = self.mng_fits.get_polyids(self.data.RA,self.data.DEC)
        self.assertTrue(numpy.all(ply == fits))
    #...

    def test_fits_areas(self):
        """Compare reading data from a .ply file with a .fits file: areas."""
        ply = self.mng.get_areas(self.data.RA,self.data.DEC)
        fits = self.mng_fits.get_areas(self.data.RA,self.data.DEC)
        self.assertTrue(numpy.all(ply == fits))
    #...

    def test_fits_weights(self):
        """Compare reading data from a .ply file with a .fits file: weights."""
        ply = self.mng.get_weights(self.data.RA,self.data.DEC)
        fits = self.mng_fits.get_weights(self.data.RA,self.data.DEC)
        self.assertTrue(numpy.all(ply == fits))
    #...

    def test_polyid(self):
        """Compare the new mangle.py with Martin White's orignal"""
        ids_new = []
        ids_orig = []
        start = time.clock()
        for x in self.data:
            ids_orig.append(self.mng_orig.polyid(x['RA'],x['DEC']))
        elapsed1 = (time.clock() - start)

        start = time.clock()
        for x in self.data:
            ids_new.append(self.mng.polyid(x['RA'],x['DEC']))
        elapsed2 = (time.clock() - start)

        self.assertEqual(ids_new,ids_orig)

        start = time.clock()
        ids_vec = self.mng.get_polyids(self.data.RA,self.data.DEC)
        elapsed3 = (time.clock() - start)
        # the results of the vector version are a numpy array, not a list
        ids_orig = numpy.array(ids_orig)
        self.assertTrue(numpy.all(ids_vec == ids_orig))
        print 'Elapsed time for polyid (orig,new):',
        print elapsed1,elapsed2
    #...

    def test_polyid_fast(self):
        """Compare the new mangle.py fast code with the C/Fortran version."""
        start = time.clock()
        ids_vec = self.mng.get_polyids(self.data.RA,self.data.DEC)
        elapsed1 = (time.clock() - start)
        ids_cmd,elapsed2 = self.do_polyid_cmd()
        self.assertTrue(numpy.all(ids_vec == ids_cmd))
        print 'Elapsed time for polyid (new_vector,C/Fortran):',
        print elapsed1,elapsed2
    #...

    def test_area(self):
        # Compare the new mangle.py with Martin White's orignal
        area_new = []
        area_orig = []
        start = time.clock()
        for x in self.data:
            area_orig.append(self.mng_orig.area(x['RA'],x['DEC']))
        elapsed1 = (time.clock() - start)

        start = time.clock()
        for x in self.data:
            area_new.append(self.mng.area(x['RA'],x['DEC']))
        elapsed2 = (time.clock() - start)

        self.assertEqual(area_new,area_orig)

        start = time.clock()
        area_new = self.mng.get_areas(self.data.RA,self.data.DEC)
        elapsed3 = (time.clock() - start)
        # the results of the vector version are a numpy array, not a list
        area_orig = numpy.array(area_orig)
        self.assertTrue(numpy.all(area_new == area_orig))
        print 'Elapsed time for area (orig,new,new_vector):',
        print elapsed1,elapsed2,elapsed3
    #...

    def test_weights(self):
        # Compare the new mangle.py with Martin White's orignal
        weight_new = []
        weight_orig = []
        start = time.clock()
        for x in self.data:
            weight_orig.append(self.mng_orig.weight(x['RA'],x['DEC']))
        elapsed1 = (time.clock() - start)

        start = time.clock()
        for x in self.data:
            weight_new.append(self.mng.weight(x['RA'],x['DEC']))
        elapsed2 = (time.clock() - start)

        self.assertEqual(weight_new,weight_orig)

        start = time.clock()
        weight_new = self.mng.get_weights(self.data.RA,self.data.DEC)
        elapsed3 = (time.clock() - start)
        # the results of the vector version are a numpy array, not a list
        weight_orig = numpy.array(weight_orig)
        self.assertTrue(numpy.all(weight_new == weight_orig))
        print 'Elapsed time for weight (orig,new,new_vector):',
        print elapsed1,elapsed2,elapsed3
    #...

    def test_getitem(self):
        """TBD: figure out a way to test these various slices."""
        mng2 = self.mng[1]
        mng2 = self.mng[1:10]
        mng2 = self.mng[[1,2,3,4,5,6,7,8,9]]
        test = self.mng.polyids < 10
        mng2 = self.mng[test]
    #...

    def test_write(self):
        """Test writing the polygon file and re-reading it."""
        self.mng.writeply(self.writetestfile)
        mng_temp = mangle.Mangle(self.writetestfile)
        numpy.testing.assert_allclose(mng_temp.polyids,self.mng.polyids)
        numpy.testing.assert_allclose(mng_temp.areas,self.mng.areas)
        numpy.testing.assert_allclose(mng_temp.weights,self.mng.weights)
        for poly1,poly2 in zip(mng_temp.polylist,self.mng.polylist):
            numpy.testing.assert_allclose(poly1,poly2)
        if 'pixels' in dir(self.mng):
            numpy.testing.assert_allclose(mng_temp.pixels,self.mng.pixels)
    #...
#...

if __name__ == '__main__':
    unittest.main()
