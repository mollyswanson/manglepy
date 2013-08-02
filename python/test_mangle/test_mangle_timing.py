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
import mangleTester

class TestMangleTiming(unittest.TestCase,mangleTester.MangleTester):
    def __init__(self, *args, **kwargs):
        """
        Set Npoints=100000 to make the test run longer to better compare timings.
        """
        mangleTester.MangleTester.__init__(self,*args,**kwargs)
        unittest.TestCase.__init__(self, *args, **kwargs)
    #...

    def setUp(self):
        """
        Set up for an individual test.

        We want to create a new mangle instance for every test,
        incase we mess it up while testing. This makes the tests take longer,
        but that's ok as it isn't part of the timed portions.
        """
        mangleTester.MangleTester.setUp(self)
        # Write out a file containing RA DEC for the commandline mangle
        outfile = file(self.mangletestfile,'w')
        outcsv = csv.writer(outfile,delimiter=' ')
        outfile.write('#RA Dec\n')
        outfile.write('#This file written by test_mangle for unittesting.\n')
        outfile.write('#It may be safely deleted.\n')
        outcsv.writerows(zip(self.data.RA,self.data.DEC))
        outfile.close()
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
#...

if __name__ == '__main__':
    unittest.main()
