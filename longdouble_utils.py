#longdouble_utils.py
"""
Basic utility functions for converting between strings and np.float128 floating point numbers.

np.float128's are actually C-style longdoubles in an 80-bit format, not full quad-precision,
(see http://en.wikipedia.org/wiki/Extended_precision) but are stored in 128 bits of memory.

Python's native methods for converting between strings and floating point numbers do
not work for np.float128's because the numbers get converted to np.float64's in the process.

Example:
import numpy as np
np.longdouble(0.1234567890123456789)

produces an output of
0.12345678901234567737

The last few digits are lost here because the string of numbers is first converted
to a double (np.float64), which can only store ~15 significant decimal digits.

To avoid this issue, this package uses the arbitrary floating point precision mpmath
http://code.google.com/p/mpmath/
to convert the string and turn it into the correct longdouble number:

import longdouble_utils.py as ld
ld.string2longdouble('0.1234567890123456789')

produces an output of
0.1234567890123456789

Here, ~19 decimal digits of precision are accurately retained as a longdouble.

This package uses an intermediary format of a 'doubledouble', which is a tuple of 2
doubles where the first double has the same value as the result of simply casting the input string
to a double, and the 2nd double contains a correction factor such that upconverting both of the
doubles to longdoubles and adding them together yields the correct longdouble value.

This technique is based on algorithm from
http://web.mit.edu/tabbott/Public/quaddouble-debian/qd-2.3.4-old/docs/qd.pdf

Functions included interconvert between strings, doubledoubles, and longdoubles:
string2doubledouble
doubledouble2string
doubledouble2longdouble
longdouble2doubledouble
string2longdouble
longdouble2string
"""
import numpy as np
import mpmath

__author__ = "Molly Swanson"
__version__ = "1.0"
__email__  = "molly.swanson@gmail.com"

def string2doubledouble(inputstr):
    """
    Convert a numerical string representing a high precision floating point number to a double-double,
    which is a tuple of two doubles that can be cast to long doubles and added together.  The first
    element in the tuple is the same as the result of just casting the string to a double, and the second
    element is a correction factor that gives the difference between the first element and the closest
    longdouble approximation to the string.  Longdoubles are assumed to be c longdoubles with 80 bit
    extended precision: http://en.wikipedia.org/wiki/Extended_precision

    Based on algorithm from http://web.mit.edu/tabbott/Public/quaddouble-debian/qd-2.3.4-old/docs/qd.pdf
    """
    mpmath.mp.prec=105
    a=mpmath.mpmathify(inputstr)
    t=mpmath.fmul((2**53+1),a)
    hi=mpmath.fsub(t,mpmath.fsub(t,a))
    lo=mpmath.fsub(a,hi,prec=11)
    lo=mpmath.fsub(a,hi)
    floathi=np.float64(hi)
    floatlo=np.float64(lo)
    doubledouble=(floathi,floatlo)
    return doubledouble     

def doubledouble2string(doubdoub,n=19):
    """
    Convert a double-double tuple into a string with n decimal digits of precision.  Default is n=19,
    which is the maximum precision that this longdouble-based double-double format can provide.
    min_fixed is the position of the leading digit such that any number with a leading
    digit less than or equal to min_fixed will be written in exponential format.
    e.g., 0.00123 has its leading digit in the -3 place, so if min_fixed>=-3
    it will be printed as 1.23e-3, and if min_fixed<=-4 it will be printed as 0.00123.
    """
    mpmath.mp.prec=105
    hi=mpmath.mpmathify(doubdoub[0])
    lo=mpmath.mpmathify(doubdoub[1])
    tot=mpmath.fadd(hi,lo)
    return mpmath.nstr(tot,n)

def doubledouble2longdouble(doubdoub):
    """
    Convert a double-double tuple into a longdouble.
    """
    longdoub=np.longdouble(doubdoub[0])+np.longdouble(doubdoub[1])
    return longdoub

def longdouble2doubledouble(longdoub):
    """
    Convert a longdouble into a double-double tuple.
    """
    floathi=np.float64(longdoub)
    floatlo=np.float64(longdoub-np.float128(floathi))
    doubledouble=(floathi,floatlo)
    return doubledouble

def string2longdouble(inputstr):
    """
    Convert a numerical string into a longdouble.  Accurately provides the full longdouble
    precision (~19 decimal digits), unlike simply doing longdouble(inputstr), which will
    internally cast inputstr to a double first.
    """
    longdoub=doubledouble2longdouble(string2doubledouble(inputstr))
    return longdoub     

def longdouble2string(longdoub,n=19):
    """
    Convert a longdouble into a numerical string with n decimal digits of precision.
    Default is n=19, which is the maximum precision that a longdouble can provide.
    min_fixed is the position of the leading digit such that any number with a leading
    digit less than or equal to min_fixed will be written in exponential format.
    e.g., 0.00123 has its leading digit in the -3 place, so if min_fixed>=-3
    it will be printed as 1.23e-3, and if min_fixed<=-4 it will be printed as 0.00123.
    """
    outputstr=doubledouble2string(longdouble2doubledouble(longdoub),n=n)
    return outputstr     
