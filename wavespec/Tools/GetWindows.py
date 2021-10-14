import numpy as np

def GetLSWindows(t,wind,slip,ngd,Ti0,Ti1,Tax=None):
	'''
	Gets the total nuimber of time windows to transform.
	
	Inputs
	======
	t: 		time array
	wind: 	window length in seconds
	slip:	amount of time to advance the window in seconds
	ngd:	number of good sections of data - from DetectGaps
	Ti0: 	starting index of every good data section - from DetectGaps
	Ti1: 	final index of every good data section - from DetectGaps
	LenW: 	set to an integer value to force a specific number of elements 
			in all windows, otherwise the window length is defined by wind
	'''

	if Tax is None:
		Tranges = t[Ti1] - t[Ti0]
		Nw = np.int32((Tranges - wind)/slip) + 2

		tax = []
		Wi0 = []
		Wi1 = []
		for i in range(0,ngd):
			t0 = t[Ti0] + np.arange(Nw[i])*slip
			t1 = t0 + wind
			tax.append(t0 + 0.5*wind)
			if Nw[i] > 0:
				i0 = np.zeros(Nw[i],dtype='int32')
				i1 = np.zeros(Nw[i],dtype='int32')
				for j in range(0,Nw[i]):
					use = np.where((t >= t0[j]) & (t < t1[j]))[0]
					if use.size == 0:
						i0[j] = -1
						i1[j] = -1
					else:
						i0[j] = use.min()
						i1[j] = use.max()
			else:
				i0 = []
				i1 = []
			Wi0.append(i0)
			Wi1.append(i1)
			
	else:

		Nw = np.array([Tax.size])
		tax = []
		Wi0 = []
		Wi1 = []		
		

		t0 = Tax - wind/2.0
		t1 = Tax + wind/2.0	
		tax.append(t0 + 0.5*wind)

		i0 = np.zeros(Nw[0],dtype='int32')
		i1 = np.zeros(Nw[0],dtype='int32')
		for j in range(0,Nw[0]):
			use = np.where((t >= t0[j]) & (t < t1[j]))[0]
			if use.size == 0:
				i0[j] = -1
				i1[j] = -1
			else:
				i0[j] = use.min()
				i1[j] = use.max()

		Wi0.append(i0)
		Wi1.append(i1)
		Nw = np.array(Nw)
		
	posWind = np.where(Nw > 0)[0]
	NwTot = np.sum(Nw[posWind]) + ngd - 1	
	
	Tax = np.zeros(NwTot,dtype='float64')
	nd = 0
	pos = 0
	for i in range(0,ngd):
		if nd > 0:
			Tax[pos] = (Tax[pos-1] + t[Ti0[i]] + wind/2.0)/2.0
			pos += 1
		
		if Nw[i] > 0:
			Tax[pos:pos+Nw[i]] = tax[i]
			pos += Nw[i]
			nd += 1		
		
			
	
	return NwTot,Nw,Wi0,Wi1,Tax

def GetFFTWindows(t,wind,slip,ngd,Ti0,Ti1,WindowUnits='s',Tax=None):
	'''
	Gets the total nuimber of time windows to transform.
	
	Inputs
	======
	t: 		time array
	wind: 	window length in seconds
	slip:	amount of time to advance the window in seconds
	ngd:	number of good sections of data - from DetectGaps
	Ti0: 	starting index of every good data section - from DetectGaps
	Ti1: 	final index of every good data section - from DetectGaps
	WindowUnits: 's'|'e' (seconds or elements)
	'''

	#get the time resolution
	Res = t[1] - t[0]

	#get the length of the window (in numbers of elements)
	#and slip
	if WindowUnits == 's':
		LenW = np.int32(np.round(wind/Res))
		LenS = np.int32(np.round(slip/Res))
	else:
		LenW = np.int32(np.round(wind))
		LenS = np.int32(np.round(slip))
	Len = Ti1 - Ti0 + 1
	
	if Tax is None:
		#get window indices
		Nw = []
		Wi0 = []
		Wi1 = []
		for i in range(0,ngd):
			nw = (Len[i] - LenW)//LenS + 1
			Nw.append(nw)
			i0 = np.arange(nw)*LenS + Ti0[i]
			i1 = i0 + LenW - 1
			Wi0.append(i0)
			Wi1.append(i1)
	else:
		Nw = np.array([Tax.size])
		tax = []
		Wi0 = []
		Wi1 = []		
		

		t0 = Tax - wind/2.0
		t1 = Tax + wind/2.0	
		tax.append(t0 + 0.5*wind)

		i0 = np.zeros(Nw[0],dtype='int32')
		i1 = np.zeros(Nw[0],dtype='int32')
		for j in range(0,Nw[0]):
			use = np.where((t >= t0[j]) & (t < t1[j]))[0]
			if use.size == 0:
				i0[j] = -1
				i1[j] = -1
			else:
				i0[j] = use.min()
				i1[j] = use.max()

		Wi0.append(i0)
		Wi1.append(i1)
		Nw = np.array(Nw)

	Nw = np.array(Nw)
	posWind = np.where(Nw > 0)[0]
	NwTot = np.sum(Nw[posWind]) + ngd - 1	
		
	
	return NwTot,LenW,LenS,Nw,Wi0,Wi1,Tax

