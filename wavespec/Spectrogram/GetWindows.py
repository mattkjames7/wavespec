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
		LenW = np.int32(np.round(wind/Res))//2

	return NwTot,LenW,Nwind

