from distutils.core import setup
from distutils.extension import Extension
#from Cython.Distutils import build_ext
from numpy.distutils.misc_util import get_numpy_include_dirs

#incdir = get_numpy_include_dirs() + ["."]
incdir = get_numpy_include_dirs()

#ext_modules = [Extension("mangle_utils", ["mangle_utils.c"], include_dirs=incdir)]
ext_modules = [Extension("mangle.mangle_utils", ["mangle/mangle_utils.c"], include_dirs=incdir)]

setup(
    name = 'mangle',
    ext_modules = ext_modules, 
    packages=['mangle']
)

