import numpy as np
from scipy.signal import convolve
from scipy.interpolate import InterpolatedUnivariateSpline

def _MakeFilter(cutoff_period,sample_freq,ftype='high'):
	
	cutoff_freq = 1.0/cutoff_period
	nyquist_freq = sample_freq/2.0
	no_nyquist = cutoff_freq/nyquist_freq
	filter_len = 3*cutoff_period*nyquist_freq
	
	if filter_len < 3:
		print('filter_len too short')
		return None
	
	#make Lanczos squared filter
	N = np.int32(filter_len)
	fltr = np.zeros(2*N-1,dtype='float64')
	fltr[N-1] = 1.0
	for i in range(0,2*N-1):
		if i != N-1:
			fltr[i] = (np.sin(np.abs(i-N+1)*np.pi/(N-1))*(N-1)/(np.pi*np.abs(i-N+1)))**2
	
	#apply cutoff factor (down 6dB at no_nyquist nyquists)
	if no_nyquist > 0:
		for i in range(0,2*N-1):
			if i != N-1:
				fltr[i] = fltr[i]*np.sin((i-N+1)*np.pi*no_nyquist)/((i-N+1)*np.pi*no_nyquist)
	
	#determine normalisation factor
	norm = (2*N-1)/np.sum(fltr)
	
	if ftype == 'high' or ftype == 'h':
		#return high pass filter
		fltr =- fltr*norm
		fltr[N-1] = fltr[N-1]+2*N-1
	else:
		#return low pass filter
		fltr = fltr*norm
	
	#normalise to length of filter
	fltr = fltr/(2*N-1)
	
	return fltr

def Filter(data,inter,high=None,low=None,KeepDC=False):
	'''
	This function performs a Lanczos-squared filter on a time series.
	
	Inputs:
		data: time series, evenly sampled.
		high: high pass cutoff period in seconds. If set equal to inter,
			  then no high pass filtering will be performed.
		low: low pass cutoff period in seconds. If set equal to inter,
			 then low pass filtering will not be performed.
		inter: time interval between time series samples in seconds.
		KeepDC: if True, the DC component of the signal will be added 
				back to the output signal.
				
	Returns:
		Filtered time series.
	'''
	if high is None:
		high = inter
	if low is None:
		low = inter
	
	
	#find bad data
	bad = np.where(np.logical_not(np.isfinite(data)))[0]
	nb = np.size(bad)
	tmpdata = np.array(data)
	l = np.size(data)


	#interpolate crap
	if nb > 0 and l-nb > 3:
		gd = np.where(np.isfinite(data))[0]
		tmp = np.arange(l)
		f = InterpolatedUnivariateSpline(tmp[gd],data[gd])
		tmpdata[bad] = f(tmp[bad])
	elif nb > 0:
		tmpdata[bad] = 0.0

	#remove DC component
	mean = np.sum(tmpdata)/l
	tmpdata -= mean


	#perform low-pass filter
	if (low > inter):
		fltr = _MakeFilter(np.float(low),1.0/np.float(inter),ftype='low')
		if not fltr is None:
			tmpdata = convolve(tmpdata,fltr)

	ts = tmpdata.size
	
	if l % 2 == 1:
		tmpdata = tmpdata[np.int(ts/2)-np.int(l/2):np.int(ts/2)+np.int(l/2)+1]
	else:
		tmpdata = tmpdata[np.int(ts/2)-np.int(l/2):np.int(ts/2)+np.int(l/2)]


	#perform high-pass filter
	if (high > inter):
		fltr = _MakeFilter(np.float(high),1.0/np.float(inter),ftype='high')
		if not fltr is None:
			tmpdata = convolve(tmpdata,fltr)		

	ts = tmpdata.size
	
	if l % 2 == 1:
		tmpdata = tmpdata[np.int(ts/2)-np.int(l/2):np.int(ts/2)+np.int(l/2)+1]
	else:
		tmpdata = tmpdata[np.int(ts/2)-np.int(l/2):np.int(ts/2)+np.int(l/2)]


	#add bad data back in
	if nb > 0:
		tmpdata[bad] = data[bad]
		
	#add DC component back in
	if KeepDC:
		tmpdata += mean
	return tmpdata
