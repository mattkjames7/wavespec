import numpy as np


def _GetFFTWindows(t,Res,wind,slip,Tax):
	'''
	Work out the start and end indices of each window in the 
	spectrogram.
	
	NOTE: while the units in thie function are expressed in s and Hz, 
	they can be substituted with m and m^-1, respetively, if doing a
	spatial FFT.
	
	Inputs
	======
	t : float
		Time array in seconds.
	Res : float
		Time resolution (s).
	wind : float
		Window length (s).
	slip : float
		Time difference from one window to the next
	Tax : None|float
		If None (default) then everything is calculated based on the 
		other input parameters automatically. If an array of times is 
		provided, the routine will do it's best to use those as the 
		middle of each window. Any from outside of the limits of t will 
		be removed.
	
	Returns
	=======
	nw : int
		The number of windows.
	lw : int
		The length of a window in elements.
	ls : int 
		The number of elements separating each window.
	i0 : int
		Array of start indices for each window.
	i1 : int
		Array of end indices for each window (i0 + lw).

	NOTE: to select the correct number of elements of the array for a
	specific window, e.g.:
	w = 3
	tw = t[i0[w]:i1[w]]
	
	'''
	
	#calculate the window length and slip length first
	lw = np.int32(np.round(wind/Res))
	ls = np.int32(np.round(slip/Res))
	
	#
	if Tax is None:
		#get the number of windows
		nw = np.int32((t.size - lw)//ls + 1)
		
		#get the start and end indices
		i0 = np.arange(nw)*ls
		i1 = i0 + lw
		
		#get the mid time
		tmid = np.float64(t[0] + Res*(lw/2.0  + i0.astype('float64')))
	else:
		#get a provisional number of windows
		nw = Tax.size
		
		#get estimates of the start and end bin numbers
		i0 = np.int32(np.round((Tax - wind/2.0 - t[0])/Res))
		i1 = i0 + lw
		
		#find out which are within the index limits of the array
		use = np.where((i0 >= 0) & (i1 <= t.size))[0]
		
		#update
		nw = use.size
		i0 = i0[use]
		i1 = i1[use]
		tmid = np.float64(t[0] + Res*(lw/2.0  + i0.astype('float64')))
		
	return nw,lw,ls,i0,i1,tmid
