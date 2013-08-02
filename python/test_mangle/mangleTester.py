"""
Helps setup and teardown of mangle tests.
"""
import numpy as np
import os

import mangle

class MangleTester(object):
    def __init__(self, *args, **kwargs):
        """
        Set up things for every test, so we don't have to regenerate 
        randoms every test.
        
        Set Npoints=100000 to make the test run longer to better compare timings.
        Set seed to change the numpy seed (default=100).
        """
        Npoints = kwargs.get('Npoints',1000)
        seed = kwargs.get('seed',100)
        np.random.seed(seed)
        # to work with the commandline version.
        self.mangletestfile = './mangle_test_data.dat'
        
        # TBD jkp: need to either include these in the mangle release
        # or create some sample polygon files (one pixelized, one not).
        # A smaller example polygon file would reduce the time it takes to run
        # the tests, since reading the files is a substantial fraction of the runtime.
        #self.polyfile = '../../geometry/mask_DR10v8_CMASS_North.ply'
        self.polyfile = '../../geometry/boss_geometry_2011_06_10.ply'
        self.polyfile_pix = '../../geometry/centerpost_mask.ply'
        
        # generate some random on-sky points
        z = np.random.uniform(-1,1,size=Npoints)
        RA = np.random.uniform(0,360,size=Npoints)
        Dec = np.degrees(np.arcsin(z))
        formats = ['f8','f8']
        names = ['RA','DEC']
        dtype = zip(names,formats)
        self.data = np.array(zip(RA,Dec),dtype=dtype).view(np.recarray)
        self.writetestfile = 'mangle_test_polys.ply'

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
        
        # TBD: just require that the fits file exist for now.
        # Should eventually have a test that compares .ply and .fits versions.
        self.fitsfile = self.polyfile.split('.ply')[0]+'.fits'
        if os.path.exists(self.fitsfile):
            self.mng_fits = mangle.Mangle(self.fitsfile)
        if os.path.exists(self.fitsfile+'.gz'):
            self.mng_fits = mangle.Mangle(self.fitsfile+'.gz')
        else:
            self.mng_fits = None
    
    def tearDown(self):
        if os.path.exists(self.mangletestfile):
            os.remove(self.mangletestfile)
        if os.path.exists(self.writetestfile):
            os.remove(self.writetestfile)
#...
