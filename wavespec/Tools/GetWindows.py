import numpy as np

def GetWindows(t,wind,slip,ngd,Ti0,Ti1,LenW=None):
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

	Res = t[1] - t[0]
	Tranges = t[Ti1] - t[Ti0]
	Nwind = np.int32((Tranges - wind)/slip) + 2
	
	posWind = np.where(Nwind > 0)[0]
	NwTot = np.sum(Nwind[posWind]) + ngd - 1	
	
	if LenW is None:
		LenW = np.int32(np.round(wind/Res))

	return NwTot,LenW,Nwind

def GetFFTWindows(t,wind,slip,ngd,Ti0,Ti1,WindowUnits='s'):
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

	Nw = np.array(Nw)
	posWind = np.where(Nw > 0)[0]
	NwTot = np.sum(Nw[posWind]) + ngd - 1	
		
	
	return NwTot,LenW,LenS,Nw,Wi0,Wi1

