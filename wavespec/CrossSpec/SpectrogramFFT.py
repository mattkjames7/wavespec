import numpy as np
from .CrossSpec import CrossSpec
from ..Tools.GetWindows import GetFFTWindows
from ..Tools.DetectGaps import DetectGaps
from ..Tools.PolyDetrend import PolyDetrend
from ..Tools.RemoveStep import RemoveStep

def SpectrogramFFT(t,v,wind,slip,WindowFunction=None,Param=None,
				FreqLim=None,Detrend=True,FindGaps=True,GoodData=None,
				Quiet=True,Threshold=0.0,WindowUnits='s',
				OneSided=True,Steps=None):
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
	Tax :	(LS only)
			An array of times at the centre of each bin.
									
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

	#find out the length of the array and 
	Tlen = np.size(t)
	if Tlen <= 1:
		return (0,0,0,0)

	#work out the time resolution - let's hope that it is constant
	Res = t[1] - t[0]

	#detect and gaps in the input data
	if FindGaps:
		ngd,Ti0,Ti1 = DetectGaps(v,GoodData)
	else:
		ngd = 1
		Ti0 = np.array([0])
		Ti1 = np.array([Tlen-1])

	#find the number of windows
	NwTot,LenW,LenS,Nw,Wi0,Wi1 = GetFFTWindows(t,wind,slip,ngd,Ti0,Ti1,WindowUnits)
	
	#find the number of frequencies
	Freq = np.arange(LenW+1,dtype='float32')/(LenW*Res)
	if OneSided:
		Freq = Freq[:LenW//2 + 1]	
	if FreqLim is None:
		Nf = Freq.size - 1		
		find = np.arange(Nf).astype('int32')
	else:
		find = np.where((Freq >= FreqLim[0]) & (Freq <= FreqLim[1]))[0]
		Nf = usef.size


	#create the output arrays
	dtype = [	('Tspec','float64'),		#mid point in time of the current window
				('Pow','float32',(Nf,)),	#Power spectra
				('Pha','float32',(Nf,)),	#phase spectra
				('Amp','float32',(Nf,)),	#Amplitude
				('Comp','complex64',(Nf,)),	#Real and Imaginary components of spectra
				('Size','int32'),			#Number of valid (finite) values used to create spectrum
				('Good','float32'),			#Fraction of good data
				('Var','float32'),]			#Variance
	out = np.recarray(Nw,dtype=dtype)
	out.fill(np.nan)
	out.nV = 0.0
	
	#Populate the time axis
	nd = 0
	pos = 0
	for i in range(0,ngd):
		if nd > 0:
			#this bit adds a load of NaNs in a gap in the middle of two good sections
			out.Tspec[pos] = (out.Tspec[pos-1] + t[Ti0[i]] + wind/2.0)/2.0
			pos+=1
		
		if Nw[i] > 0:
			out.Tspec[pos:pos+Nw[i]] = np.arange(Nw[i],dtype='float64')*slip + wind/2.0 + t[Ti0[i]]
			pos += Nw[i]
			nd += 1
				
	#loop through each good secion of the time series and FFT
	nd=0
	pos=0
	ind0 = np.arange(LenW).astype('int32')
	for i in range(0,ngd):
		if nd > 0:
			#this bit adds a load of NaNs in a gap in the middle of two good sections
			pos+=1
		
		if Nw[i] > 0:
			#loop through each window
			for j in range(0,Nw[i]):
				#get the data for this window
				ind = Wi0[i][j] + ind0
				tw = t[ind]
				vw0 = v[0][ind]
				vw1 = v[1][ind]
				
				bad = ((np.isfinite(vw0) == False) | (np.isfinite(vw1) == False)).any()
				
				#assuming everything is good, go ahead with the FFT
				if not bad:
					#remove steps and					
					#detrend if necessary
					if not Steps is None:
						vw0 = RemoveStep(tw,vw0,Steps[ind],2,5)						
						vw1 = RemoveStep(tw,vw1,Steps[ind],2,5)						
					if Detrend:
						vw0 = PolyDetrend(tw,vw0,np.int(Detrend))
						vw1 = PolyDetrend(tw,vw1,np.int(Detrend))


					power,amp,phase,fr,fi,freq = CrossSpec(tw,vw0,vw1,Freq,'FFT',WindowFunction,Param,Threshold=Threshold,Fudge=Fudge,OneSided=OneSided)
					out.Var[j+pos] = np.var(vw0) + np.var(vw1)
				
					out.Pow[j+pos] = power[find]
					out.Pha[j+pos] = phase[find]
					out.Amp[j+pos] = amp[find]
					out.Comp[j+pos] = fr[find] + 1j*fi[find]
					out.Size[j+pos] = ind.size
					out.Good[j+pos] = 1.0
				else:
					out[j+pos].Size = 0
				if not Quiet:
					print('\r{:6.2f}%'.format(100.0*np.float32(pos+j+1)/NwTot),end='')
	
			pos += Nw[i]
			nd += 1
	if not Quiet:
		print('')
			
	return NwTot,Freq,out
	
