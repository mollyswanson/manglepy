#!/usr/bin/env python
"""
Test class for the mangle version of python.

The goal is to test both the speed and correctness of improvements to
mangle.py, compared to the known-good, fast C/Fortran version.

Run with:
    python test_mangle.py
"""

# Initial Version: 2010.07.14 John K. Parejko
# 2013.06.04 jkp: updated to test pixelized masks, and do everything vs. cmdline version

import unittest
import numpy as np
import os, os.path
import pyfits
import csv
import subprocess
import time
import tempfile

import mangle

class TestMangle(unittest.TestCase):
    def __init__(self,*args, **kwargs):
        """Set up things for every test,
        so we don't have to regenerate randoms every test"""
        # TBD jkp: need to either include these in the mangle release
        # or create some sample polygon files (one pixelized, one not).
        # A smaller example polygon file would reduce the time it takes to run
        # the tests, since reading the files is a substantial fraction of the runtime.
        self.polyfile = '../../geometry/boss_geometry_2012_06_18.ply'
        self.polyfile_pix = '../../geometry/centerpost_mask.ply'
        np.random.seed(100)

        # generate some random on-sky points
        # NOTE: increase the number of points to make test run longer,
        # to better measure the commandline vs. python versions.
        Npoints = 10000
        z = np.random.uniform(-1,1,size=Npoints)
        RA = np.random.uniform(0,360,size=Npoints)
        Dec = np.degrees(np.arcsin(z))
        formats = ['f8','f8']
        names = ['RA','DEC']
        dtype = zip(names,formats)
        self.data = np.array(zip(RA,Dec),dtype=dtype).view(np.recarray)

        self.writetestfile = 'mangle_test_polys.ply'

        unittest.TestCase.__init__(self, *args, **kwargs)
    #...

    def setUp(self):
        """
        Set up for an individual test.

        We want to create a new mangle instance for every test,
        incase we mess it up while testing. This makes the tests take longer,
        but that's ok as it isn't part of the timed portions.
        """
        # set up the mangle polygon geometries
        self.mng = mangle.Mangle(self.polyfile)
        self.mng_pix = mangle.Mangle(self.polyfile_pix)
        #print 'Using',self.polyfile,'as the polygon file.'
        #print 'Using',self.polyfile,'as the pixelized polygon file.'

        # just require that the fits file exist for now.
        self.fitsfile = self.polyfile.split('.ply')[0]+'.fits'
        if os.path.exists(self.fitsfile):
            self.mng_fits = mangle.Mangle(self.fitsfile)
        else:
            self.mng_fits = None

        # Write out a file containing RA DEC for the commandline mangle
        self.mangletestfile = './mangle_test_data.dat'
        outfile = file(self.mangletestfile,'w')
        outcsv = csv.writer(outfile,delimiter=' ')
        outfile.write('#RA Dec\n')
        outfile.write('#This file written by test_mangle for unittesting.\n')
        outfile.write('#It may be safely deleted.\n')
        outcsv.writerows(zip(self.data.RA,self.data.DEC))
        outfile.close()
    #...

    def tearDown(self):
        if os.path.exists(self.mangletestfile):
            os.remove(self.mangletestfile)
        if os.path.exists(self.writetestfile):
            os.remove(self.writetestfile)
    #...

    def do_polyid_cmd(self,pix=False):
        """
        Run the commandline version of polyid.
        Return the polygons and the total time the command took.
        Set pix to use the pixelized polygon mask.
        """
        polyidoutfile = './mangle_test_polyid.dat'
        if pix:
            polyid_call = 'polyid -q '+' '.join((self.polyfile_pix,self.mangletestfile,polyidoutfile))
        else:
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
        return np.array(polyids),elapsed
    #...

    def test_fits_polyid(self):
        """Compare reading data from a .ply file with a .fits file: polyids"""
        ply = self.mng.get_polyids(self.data.RA,self.data.DEC)
        fits = self.mng_fits.get_polyids(self.data.RA,self.data.DEC)
        self.assertTrue(np.all(ply == fits))
    #...

    def test_fits_areas(self):
        """Compare reading data from a .ply file with a .fits file: areas."""
        ply = self.mng.get_areas(self.data.RA,self.data.DEC)
        fits = self.mng_fits.get_areas(self.data.RA,self.data.DEC)
        self.assertTrue(np.all(ply == fits))
    #...

    def test_fits_weights(self):
        """Compare reading data from a .ply file with a .fits file: weights."""
        ply = self.mng.get_weights(self.data.RA,self.data.DEC)
        fits = self.mng_fits.get_weights(self.data.RA,self.data.DEC)
        self.assertTrue(np.all(ply == fits))
    #...

    def test_polyid_fast(self):
        """Compare the new mangle.py fast code with the C/Fortran version."""
        start = time.clock()
        ids_vec = self.mng.get_polyids(self.data.RA,self.data.DEC)
        elapsed1 = (time.clock() - start)
        ids_cmd,elapsed2 = self.do_polyid_cmd()
        self.assertTrue(np.all(ids_vec == ids_cmd))
        print 'Elapsed time for polyid (new_vector,C/Fortran):',
        print elapsed1,elapsed2
    #...

    def test_areas(self):
        """Compare areas from mangle.py vs. areas from polyids of commandline version."""
        return
    #...
    
    def test_weights(self):
        """Compare weights from mangle.py vs. weights from polyids of commandline version."""
        return
    #...
    
    def test_pixelization(self):
        """Test how mangle.py handles a pixelized mask compared with the commandline version."""
        start = time.clock()
        ids_vec = self.mng_pix.get_polyids(self.data.RA,self.data.DEC)
        elapsed1 = (time.clock() - start)
        ids_cmd,elapsed2 = self.do_polyid_cmd(pix=True)
        self.assertTrue(np.all(ids_vec == ids_cmd))
        print 'Elapsed time for pixelized polyid (new_vector,C/Fortran):',
        print elapsed1,elapsed2
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
        np.testing.assert_allclose(mng_temp.polyids,self.mng.polyids)
        np.testing.assert_allclose(mng_temp.areas,self.mng.areas)
        np.testing.assert_allclose(mng_temp.weights,self.mng.weights)
        for poly1,poly2 in zip(mng_temp.polylist,self.mng.polylist):
            np.testing.assert_allclose(poly1,poly2)
        if 'pixels' in dir(self.mng):
            np.testing.assert_allclose(mng_temp.pixels,self.mng.pixels)
    #...
#...

if __name__ == '__main__':
    unittest.main()
