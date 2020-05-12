import numpy as np
from ..Fourier.FFT import FFT
from scipy.signal import detrend
from .GetWindows import GetWindows
from ..LombScargle.LombScargle import LombScargle
from .DetectGaps import DetectGaps
from ..CrossPhase.CrossPhase import CrossPhase

def Spectrogram(t,v,wind,slip,Freq=None,Method='FFT',WindowFunction=None,
				Param=None,Detrend=True,FindGaps=True,GoodData=None,
				Quiet=True,LenW=None,Threshold=0.0,Fudge=False,OneSided=True):
	'''
	Creates a spectogram using a sliding window.
	
	Inputs
	======
		t : time array in seconds
		v : array of values the same length as t. If using crossphase,
			this should be a list or tuple containing two arrays.
	 wind : sliding window length in seconds
	 slip : difference in time between the start of one window and the 
			next - when slip < wind, each window will have an overlap,
			when slip > wind, there will be gaps where some data will be 
			unused and when slip == wind, each window is adjacent in time.
	 Freq : a list of frequencies (Hz) to solve for - only does anything
			when using L-S
	Method : Currently either 'FFT' or 'LS'
	WindowFunction : Select a window function to apply to the data before 
			the transform, the options are: 'none','cosine-bell','hamming',
			'triangle','welch','blackman','nuttall','blackman-nuttall',
			'flat-top','cosine','gaussian'			
	Param : This parameter is used to alter some of the window functions
			(see WindowFunctions.py).
	Detrend : This will linearly detrend each time window before the 
			transform.
	FindGaps : This tells the routine to scan for data gaps, when set
			to False - the data are assumed to be perfect.
	GoodData : This can be set to a boolean array which tells the DetectGaps
			function which data points are good (True) or bad (False),
			if set to None, then any non-finite data is assumed to be bad.
	Quiet : When set to True, the function produces no stdout output; when
			False, stdout shows the progress.
	LenW : This can be set to an integer value in order to force a specific
			window length (the number of elements, as opposed to the length
			in time defined using wind)
	Threshold:	If set to a value above 0, then all values which 
			correspond to frequencies where the amplitude is less than
			Threshold are set to 0, effectively removing noise from the
			spectra.
	Fudge:	(LS Only!)
			This applies a fudge for when f == Nyquist frequency, because
			small floating point numbers have relatively large errors.
			This should only be needed if intending to reproduce a
			two-sided FFT (also, if doing this then divide A by 2 and P 
			by 4).
	OneSided: (FFT Only!)
			This should be set to remove the negative frequencies in
			the second half of the spectra. In doing so, the amplitudes
			are doubled and the powers are quadrupled.
									
	Returns
	=======
	Nw : Total number of time windows in the output array
	LenW : Length of a time window (number of elements)
	Freq : Array of frequencies in Hz.
	numpy.recarray : 
			Stores the output of the transform under the following fields:
				Tspec : Time in seconds of the middle of each time window
				Pow : Power at each frequency in each window, shape (Nw,LenW)
				Pha : Phase at each frequency in each window, shape (Nw,LenW)
				Amp : Amplitude at each frequency in each window, shape (Nw,LenW)
				Real : Real component at each frequency in each window, shape (Nw,LenW)
				Imag : Imaginary component at each frequency in each window, shape (Nw,LenW)
	'''
	isLS = 'LS' in Method
	isCP = 'CP' in Method


	#check that the frequencies exist if we are using LS
	#if Freq is None and isLS:


	#find out the length of the array and 
	Tlen = np.size(t)
	if Tlen <= 1:
		return (0,0,0,0,0,0,0,0)

	if isLS:
		#we need frequencies here, so we will assume that the data are
		#evenly spaced and that we can use the FFT frequencies
		dt,ct = np.unique((t[1:]-t[:-1]),return_counts=True)
		Res = dt[ct.argmax()]		
	else:
		Res = t[1] - t[0]

	#detect and gaps in the input data
	if FindGaps:
		ngd,Ti0,Ti1 = DetectGaps(v,GoodData)
	else:
		ngd = 1
		Ti0 = np.array([0])
		Ti1 = np.array([Tlen-1])

	#find the number of windows
	Nw,LenW,Nwind = GetWindows(t,wind,slip,ngd,Ti0,Ti1,LenW)
	
	#find the number of frequencies
	if Freq is None or not isLS:
		Freq = np.arange(LenW+1,dtype='float32')/(LenW*Res)
		if OneSided or isLS:
			Freq = Freq[:LenW//2 + 1]
	elif not Freq is None and isLS:
		df = Freq[-1] - Freq[-2]
		Freq = np.append(Freq,Freq[-1] + np.abs(df))
	Nf = Freq.size - 1		
	
	#create the output arrays
	dtype = [('Tspec','float32'),('Pow','float32',(Nf,)),('Pha','float32',(Nf,)),
			('Amp','float32',(Nf,)),('Real','float32',(Nf,)),('Imag','float32',(Nf))]
	out = np.recarray(Nw,dtype=dtype)


	out.fill(np.nan)
	
	#loop through each good secion of the time series and FFT/L-S
	nd=0
	pos=0
	for i in range(0,ngd):
		if nd > 0:
			#this bit adds a load of NaNs in a gap in the middle of two good sections
			out.Tspec[pos] = (out.Tspec[pos-1] + t[Ti0[i]] + wind/2.0)/2.0
			pos+=1
		
		if Nwind[i] > 0:
			#calculate the number of elements in this section and create
			#an array of the indices to use
			ng = Ti1[i]-Ti0[i]+1
			good = np.arange(ng) + Ti0[i]
			
			#copy the subarrays for time and v
			if isCP:
				Tv0 = v[0][good]
				Tv1 = v[1][good]
				nTv = Tv0.size
			else:
				Tv = v[good]
				nTv = Tv.size
			Tt = t[good]
			
			#output time array 
			Tax = np.arange(Nwind[i],dtype='float32')*slip + wind/2.0 + Tt[0]
			out.Tspec[pos:pos+Nwind[i]] = Tax
			
			#loop through each window
			for j in range(0,Nwind[i]):
				#indices for this current window
				if isLS:
					use = np.where((Tt >= (t[Ti0[i]] + slip*j)) & (Tt < (t[Ti0[i]] + slip*j + wind)))
				else:
					use0 = np.int32(j*slip/Res)
					use = use0 + np.arange(LenW)
								
				#this shouldn't really happen, but if the length of the array
				#doesn't match the indices, or there are dodgy values
				if np.max(use) >= nTv:
					bad = True
				else:
					if isCP:
						bad = ((np.isfinite(Tv0) == False) | (np.isfinite(Tv1) == False)).any()
					else:
						bad = (np.isfinite(Tv) == False).any()
				#assuming everything is good, go ahead with the FFT
				if not bad:
					#detrend if necessary
					if isCP:
						if Detrend:
							Tvu0 = detrend(Tv0[use])
							Tvu1 = detrend(Tv1[use])
						else:
							Tvu0 = Tv0[use]
							Tvu1 = Tv1[use]
					else:	
						if Detrend:
							Tvu = detrend(Tv[use])
						else:
							Tvu = Tv[use]
					
					if Method == 'FFT':
						power,amp,phase,fr,fi,freq = FFT(Tt[use],Tvu,WindowFunction,Param,Threshold=Threshold,OneSided=OneSided)
					elif Method == 'LS':
						power,amp,phase,fr,fi = LombScargle(Tt[use],Tvu,Freq,'C++',WindowFunction,Param,Threshold=Threshold,Fudge=Fudge)
					elif Method == 'CP-FFT':
						power,amp,phase,fr,fi,freq = CrossPhase(Tt[use],Tvu0,Tvu1,Freq,'FFT',WindowFunction,Param,Threshold=Threshold,Fudge=Fudge,OneSided=OneSided)
					elif Method == 'CP-LS':
						power,amp,phase,fr,fi,freq = CrossPhase(Tt[use],Tvu0,Tvu1,Freq,'LS',WindowFunction,Param,Threshold=Threshold,Fudge=Fudge,OneSided=OneSided)
					

					out.Pow[j+pos] = power[0:Nf]
					out.Pha[j+pos] = phase[0:Nf]
					out.Amp[j+pos] = amp[0:Nf]
					out.Real[j+pos] = fr[0:Nf]
					out.Imag[j+pos] = fi[0:Nf]
				if not Quiet:
					print('\r{:6.2f}%'.format(100.0*np.float32(pos+j+1)/Nw),end='')
	
			pos += Nwind[i]
			nd += 1
	if not Quiet:
		print('')
			
	return Nw,LenW,Freq,out
	
