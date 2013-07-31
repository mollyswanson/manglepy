#!/usr/bin/env python
"""
Profile mangle.py to find the bottlenecks.

Run with:
    python profile_mangle
load with:
    import pstats
    # when testing pixelized mask:
    pix = pstats.Stats('pix.profile')
    pix.strip_dirs().sort_stats('time').print_stats()
    # when testing non-pixelized mask:
    nopix = pstats.Stats('nopix.profile')
    nopix.strip_dirs().sort_stats('time').print_stats()
"""
import mangle
import numpy as np
import cProfile

def rand_points(Npoints,seed=101):
    """Generate Npoints randomly distributed on the sphere, in RA/Dec."""
    np.random.seed(seed)
    z = np.random.uniform(-1,1,size=Npoints)
    RA = np.random.uniform(0,360,size=Npoints)
    Dec = np.degrees(np.arcsin(z))
    formats = ['f8','f8']
    names = ['RA','DEC']
    dtype = zip(names,formats)
    return np.array(zip(RA,Dec),dtype=dtype).view(np.recarray)
#...

def prof_mangle_pix(points):
    """Profile pixelized mangle."""
    #mng = mangle.Mangle('../../geometry/centerpost_mask.ply')
    mng = mangle.Mangle('../../geometry/badfield_mask_unphot-ugriz_pix.ply')
    print mng.pixelization
    prof = cProfile.Profile()
    prof.runcall(mng.get_polyids,points.RA,points.DEC)
    return prof
#...

def prof_mangle(points):
    mng = mangle.Mangle('../../geometry/boss_geometry_2012_06_18.ply')
    prof = cProfile.Profile()
    prof.runcall(mng.get_polyids,points.RA,points.DEC)
    return prof
#...

points = rand_points(100000)

prof = prof_mangle_pix(points)
prof.print_stats('time',25)
prof.dump_stats('pix.profile')

#prof = prof_mangle(points)
#prof.print_stats('time',25)
#prof.dump_stats('nopix.profile')
