import numpy as np
from ._CFunctions import _CLombScargle,_CLombScargleSam
from ..Tools.WindowFunctions import ApplyWindowFunction

def _PyLombScargle(t,x,f,Threshold=0.0,Fudge=True):
	'''
	Python code to calculate the L-S periodogram if for whatever
	reason the C version isn't working.
	
	Inputs
	======
	t : float
		time array in seconds, length n
	x : float
		data array, length n
	f : float 
		array of frequencies to determine periodogram at, length nf
	Threshold : float
		If set to a value above 0, then all values which 
		correspond to frequencies where the amplitude is less than
		Threshold are set to 0, effectively removing noise from the
		spectra.
	Fudge : bool
		This applies a fudge for when f == Nyquist frequency, because
		small floating point numbers have relatively large errors.
		This should only be needed if intending to reproduce a
		two-sided FFT (also, if doing this then divide A by 2 and P 
		by 4).
	
	Returns
	=======
	P : float
		Periodogram power
	A : float
		Amplitude
	phi : float
		Phase in radians
	a : float
		Parameter a, see Hocke 1998.
	b : float
		Parameter b, see Hocke 1998.
	
	
	'''
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
	w = np.float32(2*np.pi*f)
	
	#calculate variance (sigma**2)
	o2 = np.float32(np.sum(x**2)/(n-1))
	
	#loop through each frequency
	for i in range(0,nf):
		if f[i] == 0:
			A[i] = np.mean(x)
			P[i] = A[i]**2.0
			phi[i] = 0.0
			a[i] = A[i]
			b[i] = 0.0
		else:
			#calculate sums to get Tau 
			ss2w = np.sum(np.sin(2*w[i]*t))
			sc2w = np.sum(np.cos(2*w[i]*t))
			
			#calculate Tau
			Tau = np.arctan2(ss2w,sc2w)/(2*w[i])
			
			#calculate w(t - Tau)
			wtT = w[i]*(t - Tau)
			
			#calculate some more sums
			syc = np.sum(x*np.cos(wtT))
			sys = np.sum(x*np.sin(wtT))
			sc2 = np.sum(np.float32(np.cos(wtT)**2))
			ss2 = np.sum(np.float32(np.sin(wtT)**2))
			
			
			#get a  and b
			rt2n = np.sqrt(2/n)
			a[i] = (0.5*rt2n*syc)/(np.sqrt(sc2))
			b[i] = (0.5*rt2n*sys)/(np.sqrt(ss2))
			if Fudge:
				if (syc < 1e-10) and (sc2 < 1e-10):
					a[i] = 0.0
				if (sys < 1e-10) and (ss2 < 1e-10):
					b[i] = 0.0
		
			#calculate the power
			P[i] = (a[i]**2 + b[i]**2)*4.0
			
			#calculate amplitude
			A[i] = np.sqrt(P[i])
			
			#calculate phase
			phi[i] = ((-np.arctan2(b[i],a[i]) -w[i]*Tau + np.pi) % (2.0*np.pi)) - np.pi
			
			c = A[i]*np.exp(1j*phi[i])
			a[i] = c.real/2.0
			b[i] = c.imag/2.0
			
	if Threshold > 0.0:
		bad = np.where(A < Threshold)[0]
		if bad.size > 0:
			A[bad] = 0.0
			P[bad] = 0.0
			phi[bad] = 0.0
			a[bad] = 0.0
			b[bad] = 0.0
	
	return P,A,phi,a,b


def LombScargle(t,x0,f,Backend='C++',WindowFunction=None,Param=None,Threshold=0.0,Fudge=False):
	'''
	Calcualtes the Lomb-Scargle periodogram for an irregularly spaced 
	time series.
	
	Inputs
	======
	t : float
		time array in seconds, length n
	x : float
		data array, length n
	f : float 
		array of frequencies to determine periodogram at, length nf
	Backend : str
		'C++' or 'Python' or 'Sam'
	WindowFunction : str
		Select a window function to apply to the data before 
		the transform, the options are: 'none','cosine-bell','hamming',
		'triangle','welch','blackman','nuttall','blackman-nuttall',
		'flat-top','cosine','gaussian'
	Param : float
		This parameter is used to alter some of the window functions
		(see WindowFunctions.py).
	Threshold : float
		If set to a value above 0, then all values which 
		correspond to frequencies where the amplitude is less than
		Threshold are set to 0, effectively removing noise from the
		spectra.
	Fudge : bool
		This applies a fudge for when f == Nyquist frequency, because
		small floating point numbers have relatively large errors.
		This should only be needed if intending to reproduce a
		two-sided FFT (also, if doing this then divide A by 2 and P 
		by 4).
	
	Returns
	=======
	P : float
		Periodogram power
	A : float
		Amplitude
	phi : float
		Phase in radians
	a : float
		Parameter a, see Hocke 1998.
	b : float
		Parameter b, see Hocke 1998.
	
	'''
	
	#preformat the input variables
	t = np.array([t],dtype='float64').flatten()
	x = np.array([x0],dtype='float64').flatten()
	f = np.array([f],dtype='float64').flatten()
	x = ApplyWindowFunction(t,x,WindowFunction,Param)
	
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
		Fudge = np.bool8(Fudge)
		Threshold = np.float64(Threshold)
		
		#Call the C++ function
		_CLombScargle(t,x,n,f,nf,P,A,phi,a,b,Threshold,Fudge)

		#return result
		return P,A,phi,a,b
	elif Backend == 'Sam':
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
		_CLombScargleSam(t.astype('float64'),x.astype('float64'),n,f.astype('float64'),nf,P,A,phi,a,b)

		#return result
		return P,A,phi,a,b
	else:
		return _PyLombScargle(t,x,f,Threshold,Fudge)
	
