import numpy as np
from ._CFunctions import _CLombScargle
from ..Tools.WindowFunctions import ApplyWindowFunction

def _PyLombScargle(t,x,f):
	
	#get array sizes
	nf = np.size(f)
	n = np.size(t)
	
	#create output arrays (Power, amplitude, phase, a, b)
	P = np.zeros((nf,),dtype='float64')
	A = np.zeros((nf,),dtype='float64')
	phi = np.zeros((nf,),dtype='float64')
	a = np.zeros((nf,),dtype='float64')
	b = np.zeros((nf,),dtype='float64')

	#convert f to omega
	w = 2*np.pi*f
	
	#calculate variance (sigma**2)
	o2 = np.sum(x**2)/(n-1)
	
	#loop through each frequency
	for i in range(0,nf):
		#calculate sums to get Tau 
		ss2w = np.sum(np.sin(2*w[i]*t))
		sc2w = np.sum(np.cos(2*w[i]*t))
		
		#calculate Tau
		Tau = np.arctan2(ss2w,sc2w)/(2*w[i])
		
		#calculate w(t - Tau)
		wtT = w[i]*(t)
		
		#calculate some more sums
		syc = np.sum(x*np.cos(wtT))
		sys = np.sum(x*np.sin(wtT))
		sc2 = np.sum(np.cos(wtT)**2)
		ss2 = np.sum(np.sin(wtT)**2)
		
		#get a  and b
		rt2n = np.sqrt(2/n)
		a[i] = rt2n*syc/np.sqrt(sc2)
		b[i] = rt2n*sys/np.sqrt(ss2)
		
	#calculate the periodogram
	a2b2 = a**2 + b**2
	P[:] = (n/(4*o2))*a2b2
	
	#calculate amplitude
	A[:] = np.sqrt(a2b2)
	
	#calculate phase
	phi[:] = -np.arctan2(b,a)
	
	return P,A,phi,a,b


def LombScargle(t,x0,f,Backend='C++',WindowFunction=None,Param=None):
	'''
	Calcualtes the Lomb-Scargle periodogram for an irregularly spaced 
	time series.
	
	Inputs
	======
	t: 		time array in seconds, length n
	x: 		data array, length n
	f: 		array of frequencies to determine periodogram at, length nf
	Backend: 'C++' or 'Python'
	WindowFunction : Select a window function to apply to the data before 
			the transform, the options are: 'none','cosine-bell','hamming',
			'triangle','welch','blackman','nuttall','blackman-nuttall',
			'flat-top','cosine','gaussian'
	Param: This parameter is used to alter some of the window functions
			(see WindowFunctions.py).
	
	Returns
	=======
	P:		Periodogram power
	A: 		Amplitude
	phi:	Phase in radians
	a:		Parameter a, see Hocke 1998.
	b:		Parameter b, see Hocke 1998.
	
	'''
	
	#preformat the input variables
	t = np.array([t],dtype='float64').flatten()
	x = np.array([x0],dtype='float64').flatten()
	f = np.array([f],dtype='float64').flatten()
	x = ApplyWindowFunction(x,WindowFunction,Param)
	
	if Backend == 'C++':
		#get array sizes
		nf = np.int32(np.size(f))
		n = np.int32(np.size(t))
		
		#create output arrays (Power, amplitude, phase, a, b)
		P = np.zeros((nf,),dtype='float64')
		A = np.zeros((nf,),dtype='float64')
		phi = np.zeros((nf,),dtype='float64')
		a = np.zeros((nf,),dtype='float64')
		b = np.zeros((nf,),dtype='float64')
		
		#Call the C++ function
		_CLombScargle(t,x,n,f,nf,P,A,phi,a,b)
		
		#return result
		return P,A,phi,a,b
	else:
		return _PyLombScargle(t,x,f)
	
