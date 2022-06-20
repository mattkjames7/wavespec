import numpy as np


def _GetLSWindows(t,wind,slip,Tax):
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
	
	#get start and end times first
	T0 = t[0]
	T1 = t[-1]
	
	
	if Tax is None:
		#get the number of windows
		nw = np.int32((T1 - T0 - wind)/slip + 1)
		
		#get the mid time
		tmid = np.float64(T0 + wind/2.0 + np.arange(nw)*slip)
	else:
		nw = Tax.size
		tmid = np.float64(Tax)
		
	#get the start indices
	i0 = np.zeros(nw,dtype='int64') - 1
	i1 = np.zeros(nw,dtype='int64') - 1
	# s = 0
	# for i in range(0,nw):
		# t0 = tmid[i] - wind/2.0
		# t1 = tmid[i] + wind/2.0
		# for j in range(s,t.size):
			# if (t[j] >= t0) and (t[j] < t1):
				# s = j
				# i0[i] = j
				# break
			# if (t[j] >= t1):
				# break
	# e = t.size
	# for i in range(0,nw):
		# t0 = tmid[i] - wind/2.0
		# t1 = tmid[i] + wind/2.0
		# for j in range(e-1,-1,-1):
			# if (t[j] >= t0) and (t[j] < t1):
				# e = j+1
				# i1[i] = j+1
				# break
			# if (t[j] < t0):
				# break
	# good = np.where((i0 > -1) & (i1 > -1))[0]
	# i0 = i0[good]
	# i1 = i1[good]
	# tmid = tmid[good]
	# nw = good.size
	
	#need a faster way of doing this
	for i in range(0,nw):
		t0 = tmid[i] - wind/2.0
		t1 = tmid[i] + wind/2.0
		use = np.where((t >= t0) & (t < t1))[0]
		if use.size > 2:
			i0[i] = use[0]
			i1[i] = use[-1]+1
	good = np.where((i0 > -1) & (i1 > -1))[0]
	i0 = i0[good]
	i1 = i1[good]
	tmid = tmid[good]
	nw = good.size		
		
	return nw,i0,i1,tmid
