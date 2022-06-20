import numpy as np
from .FFT import FFT
from ._GetFFTWindows import _GetFFTWindows
from ..Tools.PolyDetrend import PolyDetrend
from ..Tools.RemoveStep import RemoveStep

def Spectrogram(t,v,wind,slip,**kwargs):
	'''
	Creates a spectogram using a sliding window.
	
	NOTE: The descriptions below use the units s and Hz for a Fourier
	transform of time-series data. For spatial FFTs replace with 
	appropriate units (e.g. m and m^-1).
	
	Inputs
	======
	t : float
		Time array in seconds - must be equally spaced!
	v : float
		Time-series data to be transformed.
	wind : float
		Sliding window length in seconds.
	slip : float
		Difference in time between the start of one window and the 
		next - when slip < wind, each window will have an overlap,
		when slip > wind, there will be gaps where some data will be 
		unused and when slip == wind, each window is adjacent in time.
		This parameter will be ignored if the Tax keyword is set.
	
	Keyword Arguments
	=================
	FreqLim : None|float
			2-element array denoting the minimum and maximum frequency 
			to use.
	WindowFunction : None|str
			Select a window function to apply to the data before the 
			transform, the options are: 
			'none','cosine-bell','hamming','triangle','welch',
			'blackman','nuttall','blackman-nuttall','flat-top','cosine',
			'gaussian'			
	Param : float
			This parameter is used to alter some of the window functions
			(see WindowFunctions.py).
	Detrend : bool|int
			This will detrend using a polynomial (the degree of which is
			given by the Detrend keyword, e.g. Detrend=2 for quadratic)
			each time window before the transform.
	GoodData : bool
			This can be set to a boolean array which tells the function 
			which data points are good (True) or bad (False), if set to 
			None, then any non-finite data is assumed to be bad.
	Quiet : bool
			When set to True, the function produces no stdout output; 
			when False, stdout shows the progress.
	Threshold:	If set to a value above 0, then all values which 
			correspond to frequencies where the amplitude is less than
			Threshold are set to 0, effectively removing noise from the
			spectra.
	OneSided: 
			This should be set to remove the negative frequencies in
			the second half of the spectra. In doing so, the amplitudes
			are doubled and the powers are quadrupled.
	Tax :	float
			An array of times at the centre of each bin. This will force
			the program to use a specific set of windows.
									
	Returns
	=======
	Nw : int 
		Total number of time windows in the output array
	Freq : float
		Array of frequencies in Hz.
	out : numpy.recarray
		Stores the output of the transform under the following fields:
			Tspec : float
				Time in seconds of the middle of each time window
			Pow : float
				Power at each frequency in each window, shape (Nw,LenW)
			Pha : float
				Phase at each frequency in each window, shape (Nw,LenW)
			Amp : float 
				Amplitude at each frequency in each window, shape (Nw,LenW)
			Real : float
				Real component at each frequency in each window, shape (Nw,LenW)
			Imag : float
				Imaginary component at each frequency in each window, shape (Nw,LenW)
			Var : float
				Variance of dat within each window
			Good : float
				Fraction of good data in each window.
			Size : int
				Number of elements in window.
	'''
	
	Tax = kwargs.get('Tax',None)
	WindowFunction = kwargs.get('WindowFunction',None)
	Param = kwargs.get('Param',None)
	FreqLim = kwargs.get('FreqLim',None)
	Detrend = kwargs.get('Detrend',True)
	GoodData = kwargs.get('GoodData',None)
	Quiet = kwargs.get('Quiet',True)
	Threshold = kwargs.get('Threshold',0.0)
	OneSided = kwargs.get('OneSided',True)


	#find out the length of the array and 
	Tlen = np.size(t)
	if Tlen <= 1:
		return None

	#work out the time resolution - let's hope that it is constant
	Res = t[1] - t[0]

	#detect good/bad data
	if GoodData is None:
		good = np.isfinite(v)
	else:
		good = GoodData
			
	#get the windows and their indices etc.
	Nw,LenW,ls,i0,i1,Tmid = _GetFFTWindows(t,Res,wind,slip,Tax=Tax)


	#find the number of frequencies
	Freq = np.arange(LenW+1,dtype='float64')/(LenW*Res)
	if OneSided:
		Freq = Freq[:LenW//2 + 1]	
	if FreqLim is None:
		Nf = Freq.size - 1		
		find = np.arange(Nf).astype('int32')
	else:
		find = np.where((Freq >= FreqLim[0]) & (Freq <= FreqLim[1]))[0]
		if find[-1] == Freq.size-1:
			find = find[:-1]
		Freq = np.append(Freq[find],Freq[find[-1]+1])
		Nf = find.size


	#create the output arrays
	dtype = [	('Tspec','float64'),		#mid point in time of the current window
				('Pow','float64',(Nf,)),	#Power spectra
				('Pha','float64',(Nf,)),	#phase spectra
				('Amp','float64',(Nf,)),	#Amplitude
				('Comp','complex128',(Nf,)),	#Real and Imaginary components of spectra
				('Size','int64'),			#Number of valid (finite) values used to create spectrum
				('Good','float64'),			#Fraction of good data
				('Var','float64'),]			#Variance
	out = np.recarray(Nw,dtype=dtype)
	out.fill(np.nan)
	out.nV = 0.0
	out.Tspec = Tmid

	#loop through each window and FFT
	ind0 = np.arange(LenW).astype('int32')
	for i in range(0,Nw):
		#get the data for this window
		ind = np.arange(i0[i],i1[i])
		tw = t[ind]
		vw = v[ind]
		gd = good[ind]
		ngd = np.sum(gd)	
		bad = ngd != LenW

		#assuming everything is good, go ahead with the FFT
		if not bad:
			#detrend if necessary
	
			if Detrend:
				vw = PolyDetrend(tw,vw,np.int(Detrend))

					
			power,amp,phase,fr,fi,freq = FFT(tw,vw,WindowFunction,Param,
								Threshold=Threshold,OneSided=OneSided)
			out.Var[i] = np.var(vw)
			out.Pow[i] = power[find]
			out.Pha[i] = phase[find]
			out.Amp[i] = amp[find]
			out.Comp[i] = fr[find] + 1j*fi[find]
			out.Size[i] = LenW
			out.Good[i] = 1.0
		else:
			out[i].Var = np.nanvar(vw)
			out[i].Size = LenW
			out[i].Good = ngd/LenW
		if not Quiet:
			print('\r{:6.2f}%'.format(100.0*np.float32(i+1)/Nw),end='')
	

	if not Quiet:
		print('')
			
	return Nw,Freq,out
	
