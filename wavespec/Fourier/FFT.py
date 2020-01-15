import numpy as np
from ..Tools.WindowFunctions import ApplyWindowFunction

def FFT(t,x,WindowFunction=None,Param=None):
	'''
	Performs a Fast Fourier Transform (FFT) on time series data.
	
	Inputs
	======
	t:		time array
	x:		time series data
	WindowFunction : Select a window function to apply to the data before 
			the transform, the options are: 'none','cosine-bell','hamming',
			'triangle','welch','blackman','nuttall','blackman-nuttall',
			'flat-top','cosine','gaussian'			
	Param : This parameter is used to alter some of the window functions
			(see WindowFunctions.py).
			
	Returns
	=======
	power: array of powers for each frequency
	phase: array of phases for each frequency
	freq: array of frequencies
	fr: real component of the FFT
	fi: imaginary component of the FFT
	'''
	
	#find out the array length (l) and the interval between each value (i)
	l = np.size(x)
	i = t[1] - t[0]
	
	#firstly apply the window function to the time series
	v = ApplyWindowFunction(x,WindowFunction,Param)
	
	#now to do the FFT
	f = np.fft.fft(v)/l
	fi = np.imag(f)
	fr = np.real(f)

	#calculate power
	power = (fi**2 + fr**2)
	
	#frequency
	freq = np.arange(l+1,dtype='float32')/(np.float32(l*i))
	
	#phase in radians
	phase = np.arctan2(fi,fr)
	
	return (power,phase,freq,fr,fi)
