import numpy as np
from ..LombScargle.LombScargle import LombScargle
from ..Fourier.FFT import FFT


def CrossSpec(t,x,y,Freq=None,Method='FFT',WindowFunction=None,
				Param=None,Threshold=0.0,Fudge=False,OneSided=False):
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
	t : float
		Time axis in seconds
	x : float
		Time series data
	y :	float
		Another time series data
	Freq : float
		List of frequencies to solve for (only needed for LS)
	Method : str
		'LS' or 'FFT'

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
	Fudge : bool
			This applies a fudge for when f == Nyquist frequency, because
			small floating point numbers have relatively large errors.
			This should only be needed if intending to reproduce a
			two-sided FFT (also, if doing this then divide A by 2 and P 
			by 4). (only for LS)						
				
	Returns
	=======
	Pxy : float
		Power
	Axy : float
		Amplitude
	phixy : float
		Phase
	Pxyr : float
		Real
	Pxyi : float
		Imaginary
	Freq : float
		Frequency
	
	'''
	#check that the frequencies exist if we are using LS
	if Freq is None and Method == 'LS':
		print('Please set the Freq keyword before using the LS method')
		return
	
	#Perform FFT/LS
	if Method == 'FFT':
		_,_,_,Xr,Xi,Freq = FFT(t,x,WindowFunction,Param,Threshold=Threshold,OneSided=OneSided)
		_,_,_,Yr,Yi,Freq = FFT(t,y,WindowFunction,Param,Threshold=Threshold,OneSided=OneSided)
	elif Method == 'LS':
		_,_,_,Xr,Xi = LombScargle(t,x,Freq,'C++',WindowFunction,Param,Threshold=Threshold,Fudge=Fudge)
		_,_,_,Yr,Yi = LombScargle(t,y,Freq,'C++',WindowFunction,Param,Threshold=Threshold,Fudge=Fudge)
	
	
	
	#calculate crossphase power
	Pxyr = Xr*Yr
	Pxyi = -Xi*Yi
	Pxy = np.abs(Pxyi**2 + Pxyr**2)
	
	#calculate phase
	phixy = np.arctan2(Pxyi,Pxyr)	
	
	#calculate amplitude
	Axy = np.sqrt(Pxy)
	
	return Pxy,Axy,phixy,Pxyr,Pxyi,Freq
