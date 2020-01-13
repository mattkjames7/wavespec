import numpy as np
from .WindowFunctions import ApplyWindowFunction

def FFT(t,x,WindowFunction=None,Params=None):

	#find out the array length (l) and the interval between each value (i)
	l = np.size(v)
	i = t[1] - t[0]
	
	#firstly apply the window function to the time series
	v = ApplyWindowFunction(x,WindowFunction,Params)
	
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
