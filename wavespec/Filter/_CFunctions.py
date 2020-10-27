from .. import Globals
import ctypes as ct
import numpy as np

#define some dtypes
c_bool = ct.c_bool
c_int = ct.c_int
c_float = ct.c_float
c_double = ct.c_double
c_int_ptr = np.ctypeslib.ndpointer(ct.c_int,flags="C_CONTIGUOUS")
c_float_ptr = np.ctypeslib.ndpointer(ct.c_float,flags="C_CONTIGUOUS")
c_double_ptr = np.ctypeslib.ndpointer(ct.c_double,flags="C_CONTIGUOUS")

#check if the libfilter.so is loaded
if Globals.libF is None:
	Globals.libF = ct.CDLL(Globals.libFfile)
libF = Globals.libF


#define some of the functions
_CIFilter = libF.IFilter
_CIFilter.restype = None
_CIFilter.argtypes = [	c_int,			#n - the number of elements in the time array
						c_float_ptr,	#time array
						c_float_ptr,	#time series data
						c_float,		#cutoff frequency
						c_int,			#integer filter type
						c_float_ptr]	#output filtered data
						
_CLanczosKernelLP = libF.LanczosKernelLP
_CLanczosKernelLP.restype = None
_CLanczosKernelLP.argtypes = [	c_int,			#element at centre of window
								c_int,			#number of elements in time series
								c_float_ptr,	#time array
								c_float,		#cutoff frequency
								c_float_ptr]	#output kernel array

_CLanczosKernelHP = libF.LanczosKernelHP
_CLanczosKernelHP.restype = None
_CLanczosKernelHP.argtypes = [	c_int,			#element at centre of window
								c_int,			#number of elements in time series
								c_float_ptr,	#time array
								c_float,		#cutoff frequency
								c_float_ptr]	#output kernel array




