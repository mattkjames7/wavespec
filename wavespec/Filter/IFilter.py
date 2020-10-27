import numpy as np
from .Lanczos import Lanczos
from ._CFunctions import _CIFilter


def _IFilterPy(t,x,fc,ftype):
	
	#get the number of elements
	n = t.size
	
	#create the output array
	out = np.zeros(n,dtype=x.dtype)
	
	#loop through each one at a time
	for i in range(0,n):
		#get the kernel at time t[i]
		kernel = Lanczos(t[i],t,fc,ftype)
		
		#now convolve them
		out[i] = np.sum(kernel*x)
		
	return out

def IFilter(t,x,fc,ftype='low',Backend='C++'):
	'''
	Irregular Lanczos squared filter. I think...
	
	
	'''
	#pick out good and bad data
	good = np.isfinite(x)
	gd = np.where(good)[0]
	bd = np.where(good == False)[0]
	tg = t[gd]
	xg = x[gd]
	
	
	#if we are using the python backend
	if Backend == 'Python':
		return _IFilterPy(t,x,fc,ftype)
		
	#otherwise use C++
	_n = np.int32(tg.size)
	_t = np.array(tg).astype('float32')
	_x = np.array(xg).astype('float32')
	_fc = np.float32(fc)
	if ftype == 'high':
		_ftype = np.int32(1)
	else:
		_ftype = np.int32(-1)
	_out = np.zeros(_n,dtype='float32')
	
	#call the C++ function
	_CIFilter(_n,_t,_x,_fc,_ftype,_out)
	
	#put the bad bits back in
	out = np.zeros(t.size,dtype='float32') + np.nan
	out[gd] = _out
	
	return out	
	
	

		

