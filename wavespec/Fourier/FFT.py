import numpy as np
from ..Tools.WindowFunctions import ApplyWindowFunction

def FFT(t,x,WindowFunction=None,Param=None,Threshold=0.0,OneSided=False):
	'''
	Performs a Fast Fourier Transform (FFT) on time series data.
	
	Inputs
	======
	t :	float
		time array
	x : float
		time series data
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
	OneSided : bool 
			This should be set to remove the negative frequencies in
			the second half of the spectra. In doing so, the amplitudes
			are doubled and the powers are quadrupled.

									
	Returns
	=======
	power : float
		array of powers for each frequency
	A : float
		array of amplitudes
	phase : float
		array of phases for each frequency
	fr : float
		real component of the FFT
	fi : float
		imaginary component of the FFT
	freq : float
		array of frequencies
	'''
	
	#find out the array length (l) and the interval between each value (i)
	l = np.size(x)
	i = t[1] - t[0]
	
	#firstly apply the window function to the time series
	v = ApplyWindowFunction(t,x,WindowFunction,Param)
	
	#now to do the FFT
	f = np.fft.fft(v)/l
	fi = np.imag(f)
	fr = np.real(f)

	#frequency
	freq = np.arange(l+1,dtype='float32')/(np.float32(l*i))

	#Amplitude
	A = np.abs(f)
	
	if Threshold > 0.0:
		#set everything with A < Threshold to 0 to remove noise
		bad = np.where(A < Threshold)[0]
		A[bad] = 0.0
		fr[bad] = 0.0
		fi[bad] = 0.0

	if OneSided:
		#reduce size of the arrays
		freq = freq[:l//2 + 1]
		fr = fr[:l//2]
		fi = fi[:l//2]
		A = A[:l//2]
		
		#multiply amplitudes by 2 except for DC
		A[1:] *= 2.0
	
	#calculate power
	power = A**2
	
	#phase in radians
	phase = np.arctan2(fi,fr)
	
	return (power,A,phase,fr,fi,freq)
