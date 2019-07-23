from .. import Globals
import ctypes as ct
import numpy as np

#define some dtypes
c_int = ct.c_int
c_float = ct.c_float
c_double = ct.c_double
c_int_ptr = np.ctypeslib.ndpointer(ct.c_int,flags="C_CONTIGUOUS")
c_float_ptr = np.ctypeslib.ndpointer(ct.c_float,flags="C_CONTIGUOUS")
c_double_ptr = np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")


#check if the liblombscargle.so is loaded
if Globals.libLS is None:
	Globals.libLS = ct.CDLL(Globals.libLSfile)
libLS = Globals.libLS

#define functions used from the library
_CLombScargle = libLS.LombScargle
_CLombScargle.argtypes = [c_double_ptr, c_double_ptr, c_int, c_double_ptr, c_int, c_double_ptr, c_double_ptr, c_double_ptr, c_double_ptr, c_double_ptr]
_CLombScargle.restype = None
