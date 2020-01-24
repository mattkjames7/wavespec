import numpy as np
from ..LombScargle.LombScargle import LombScargle
from ..Fourier.FFT import FFT


def CrossPhase(t,x,y,Freq=None,Method='FFT',WindowFunction=None,Param=None):
	'''
	This procedure will perform crossphase analysis of two time series
	using the method described in Waters et al., 1991 
	(https://doi.org/10.1016/0032-0633(91)90052-C):
	
	Power
	P_xy = | X* x Y | = A_xy + B_xy i
	
	Phase
	phi = arctan(B_xy/A_xy)
	
	
	Inputs
	======
	t:	Time axis in seconds
	x:	Time series data
	y:	Another time series data
	Freq:	List of frequencies to solve for (only needed for LS)
	Method:	'LS' or 'FFT'
	WindowFunction:	Define the window function to be used on both time series
	Param:	Window function parameter.
	
	Returns
	=======
	Pxy:	Power
	Axy:	Amplitude
	phixy:	Phase
	Pxyr:	Real
	Pxyi:	Imaginary
	Freq:	Frequency
	
	'''
	#check that the frequencies exist if we are using LS
	if Freq is None and Method == 'LS':
		print('Please set the Freq keyword before using the LS method')
		return
	
	#Perform FFT/LS
	if Method == 'FFT':
		_,_,Freq,Xr,Xi = FFT(t,x,WindowFunction,Param)
		_,_,Freq,Yr,Yi = FFT(t,y,WindowFunction,Param)
	elif Method == 'LS':
		_,_,_,Xr,Xi = LombScargle(t,x,Freq,'C++',WindowFunction,Param)
		_,_,_,Yr,Yi = LombScargle(t,y,Freq,'C++',WindowFunction,Param)
	
	
	
	#calculate crossphase power
	Pxyr = Xr*Yr
	Pxyi = -Xi*Yi
	Pxy = Pxyi**2 + Pxyr**2
	
	#calculate phase
	phixy = np.arctan2(Pxyi,Pxyr)	
	
	#calculate amplitude
	Axy = np.sqrt(Pxy)
	
	return Pxy,Axy,phixy,Pxyr,Pxyi,Freq
