import numpy as np
from ..Fourier.FFT import FFT
from scipy.signal import detrend
from .GetWindows import GetWindows
from ..LombScargle.LombScargle import LombScargle
from .DetectGaps import DetectGaps

def Spectrogram(t,v,wind,slip,Freq=None,Method='FFT',WindowFunction=None,Param=None,Detrend=True,FindGaps=True,GoodData=None,Quiet=True,LenW=None):
	'''
	Creates a spectogram using a sliding window.
	
	Inputs
	======
		t : time array in seconds
		v : array of values the same length as t
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
	LenW : This can be set to an integer value in order to for a specific
			window length (the number of elements, as opposed to the length
			in time defined using wind)
			
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
	#check that the frequencies exist if we are using LS
	if Freq is None and Method == 'LS':
		print('Please set the Freq keyword before using the LS method')
		return


	#find out the length of the array and 
	Tlen = np.size(t)
	if Tlen <= 1:
		return (0,0,0,0,0,0,0,0)

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
	if Method == 'LS':
		LenW = np.size(Freq)
	
	#create the output arrays
	dtype = [('Tspec','float32'),('Pow','float32',(LenW,)),('Pha','float32',(LenW,)),
			('Amp','float32',(LenW,)),('Real','float32',(LenW,)),('Imag','float32',(LenW))]
	out = np.recarray(Nw,dtype=dtype)
	if Method != 'LS':
		Freq = ((np.arange(LenW*2+1,dtype='float32')/(LenW*2))/Res)[0:LenW+1]
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
			Tv = v[good]
			Tt = t[good]
			
			#output time array 
			Tax = np.arange(Nwind[i],dtype='float32')*slip + wind/2.0 + Tt[0]
			out.Tspec[pos:pos+Nwind[i]] = Tax
			
			#loop through each window
			for j in range(0,Nwind[i]):
				#indices for this current window
				if Method == 'FFT':
					use0 = np.int32(j*slip/Res)
					use = use0 + np.arange(LenW*2)
				else:
					use = np.where((Tt >= (t[Ti0[i]] + slip*j)) & (Tt < (t[Ti0[i]] + slip*j + wind)))
				
				#this shouldn't really happen, but if the length of the array
				#doesn't match the indices, or there are dodgy values
				if np.max(use) >= np.size(Tv):
					bad = True
				else:
					bad = (np.isfinite(Tv) == False).any()
					
				#assuming everything is good, go ahead with the FFT
				if not bad:
					#detrend if necessary
					if Detrend:
						Tvu = detrend(Tv[use])
					else:
						Tvu = Tv[use]
					
					if Method == 'FFT':
						power,phase,freq,fr,fi = FFT(Tt[use],Tvu,WindowFunction,Param)
						amp = np.sqrt(power)
					elif Method == 'LS':
						power,amp,phase,fr,fi = LombScargle(Tt[use],Tvu,Freq,'C++',WindowFunction,Param)
					
					out.Pow[j+pos] = power[0:LenW]
					out.Pha[j+pos] = phase[0:LenW]
					out.Amp[j+pos] = amp[0:LenW]
					out.Real[j+pos] = fr[0:LenW]
					out.Imag[j+pos] = fi[0:LenW]
				if not Quiet:
					print('\r{:6.2f}%'.format(100.0*np.float32(pos+j+1)/Nw),end='')
	
			pos += Nwind[i]
			nd += 1
	if not Quiet:
		print('')
			
	return Nw,LenW,Freq,out
	
