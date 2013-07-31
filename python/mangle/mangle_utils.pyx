"""
Cythonized version of mangle polygon calculations.
"""
#Initial Version: John Parejko, Yale 2011

import numpy as np
cimport numpy as np
cimport cython

cdef extern from "math.h" :
  double fabs(double)
cdef extern from "math.h" :
  long double fabsl(long double)

# original, vectorized python version
##    def incap_vec(self,cap,x0,y0,z0):
##        """
##        Returns True for each (theta,phi) that lies in the cap specified by the
##        4-vector "cap" containing (x,y,z,cm), and False for the rest.
##        """
##        cd = 1.0-cap[0]*x0-cap[1]*y0-cap[2]*z0
##        return ((cap[3] < 0.0) & (cd>fabs(cap[3]))) | ((cap[3] > 0.0) & (cd<cap[3]))

# Turn off error checking!!
@cython.boundscheck(False)
@cython.wraparound(False)
def _incap(np.ndarray[np.float64_t,ndim=1] cap,
                np.ndarray[np.float64_t,ndim=1] x,
                np.ndarray[np.float64_t,ndim=1] y,
                np.ndarray[np.float64_t,ndim=1] z):
    """Returns True for each (x,y,z) that lies in the cap specified
    by the 4-vector "cap" containing (x,y,z,cm), and False otherwise."""
    cdef int ndim
    cdef double cd
    ndim = x.shape[0]
    cdef np.ndarray[np.int8_t,ndim=1] result
    result = np.empty(ndim, dtype=np.int8)
    
    for i in range(ndim):
        cd = 1.0 - cap[0]*x[i] - cap[1]*y[i] - cap[2]*z[i]
        if cap[3] < 0.0:
            if cd > fabs(cap[3]):
                result[i] = 1
            else:
                result[i] = 0
        else:
            if cd < cap[3]:
                result[i] = 1
            else:
                result[i] = 0
        # the above is equivalent to this, but short-circuited,
        # and thus slightly faster
        #if ((cap[3] < 0.0) & (cd>fabs(cap[3]))) | ((cap[3] > 0.0) & (cd<cap[3])):
        #    result[i] = 1
        #else:
        #    result[i] = 0
    return result
#...

def _incapl(np.ndarray[np.longdouble_t,ndim=1] cap,
                np.ndarray[np.longdouble_t,ndim=1] x,
                np.ndarray[np.longdouble_t,ndim=1] y,
                np.ndarray[np.longdouble_t,ndim=1] z):
    """Returns True for each (x,y,z) that lies in the cap specified
    by the 4-vector "cap" containing (x,y,z,cm), and False otherwise."""
    cdef int ndim
    cdef long double cd
    ndim = x.shape[0]
    cdef np.ndarray[np.int8_t,ndim=1] result
    result = np.empty(ndim, dtype=np.int8)
    
    for i in range(ndim):
        cd = 1.0 - cap[0]*x[i] - cap[1]*y[i] - cap[2]*z[i]
        if cap[3] < 0.0:
            if cd > fabsl(cap[3]):
                result[i] = 1
            else:
                result[i] = 0
        else:
            if cd < cap[3]:
                result[i] = 1
            else:
                result[i] = 0
        # the above is equivalent to this, but short-circuited,
        # and thus slightly faster
        #if ((cap[3] < 0.0) & (cd>fabsl(cap[3]))) | ((cap[3] > 0.0) & (cd<cap[3])):
        #    result[i] = 1
        #else:
        #    result[i] = 0
    return result
#...

def _inpoly_vec(np.ndarray[np.longdouble_t,ndim=1] cap,
                np.ndarray[np.longdouble_t,ndim=1] x,
                np.ndarray[np.longdouble_t,ndim=1] y,
                np.ndarray[np.longdouble_t,ndim=1] z):
    return None
#...
