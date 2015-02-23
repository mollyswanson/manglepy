#!/usr/bin/env python
"""
Test class for the mangle version of python. Tests atomic operations, but 
does not compare with commandline version.

Run with:
    python test_mangle.py
or give --verbose to get a print out of each test as it is run.
"""

from __future__ import print_function
import unittest
import numpy as np
import pyfits

import mangle

import mangleTester

class TestMangle(mangleTester.MangleTester,unittest.TestCase):
    def __init__(self, *args, **kwargs):
        """Set up things for every test,
        so we don't have to regenerate randoms every test"""
        mangleTester.MangleTester.__init__(self,*args,**kwargs)
        unittest.TestCase.__init__(self, *args, **kwargs)
    
    def runTest(self):
        pass
    
    def test_slice_weight(self):
        """Test slicing on polygon weights."""
        test = self.mng.weights > 0
        result1 = self.mng.get_polyids(self.data.RA,self.data.DEC)
        # crude way to elimante "unwanted" polygons.
        for i,x in enumerate(result1):
            if x not in self.mng.polyids[test]:
                result1[i] = -1
        mng2 = self.mng[test]
        result2 = mng2.get_polyids(self.data.RA,self.data.DEC)
        # the polyids themselves will be different, so we have to test on
        # the actual polygon identities.
        test1 = (result1 > 0)
        test2 = (result2 > 0)
        # Did everthing match that should have?
        self.assertTrue(np.all(test1 == test2))
        # Do the polygons themselves match 1-1?
        for x,y in zip(self.mng.polylist[result1[test1]],mng2.polylist[result2[test2]]):
            self.assertTrue(np.all(x == y))
        
    def test_slice_generic(self):
        """Test creating a generic index slice."""
        mng2 = self.mng[1]
        mng2 = self.mng[1:10]
        mng2 = self.mng[[1,2,3,4,5,6,7,8,9]]
    #...
    
    #@unittest.skip('skipping fits tests')
    def test_fits_polyid(self):
        """Compare reading data from a .ply file with a .fits file: polyids"""
        ply = self.mng.get_polyids(self.data.RA,self.data.DEC)
        fits = self.mng_fits.get_polyids(self.data.RA,self.data.DEC)
        self.assertTrue(np.all(ply == fits))
    #...

    #@unittest.skip('skipping fits tests')
    def test_fits_areas(self):
        """Compare reading data from a .ply file with a .fits file: areas."""
        ply = self.mng.get_areas(self.data.RA,self.data.DEC)
        fits = self.mng_fits.get_areas(self.data.RA,self.data.DEC)
        self.assertTrue(np.all(ply == fits))
    #...

    #@unittest.skip('skipping fits tests')
    def test_fits_weights(self):
        """Compare reading data from a .ply file with a .fits file: weights."""
        ply = self.mng.get_weights(self.data.RA,self.data.DEC)
        fits = self.mng_fits.get_weights(self.data.RA,self.data.DEC)
        self.assertTrue(np.all(ply == fits))
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
