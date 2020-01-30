import numpy as np


def DetectWavePeaks(t,f,P,Threshold,LargestPeak=False):
	'''
	Detects the peaks in wave power from a FFT or LS spectrogram.
	'''
	
	nt = t.size
	nf = f.size
	
	#find the peaks
	ispeak = np.isfinite(P[:,1:-1]) & (P[:,1:-1] > P[:,0:-2]) & (P[:,1:-1] > P[:,2:]) & (P[:,1:-1] >= Threshold)
	
	if LargestPeak:
		#select largest peak
		pkpwr = P[:,1:-1]*np.int32(ispeak)
		mx = np.max(pkpwr,axis=1)
		for i in range(0,nt):
			ispeak[i] = pkpwr[i] == mx[i]
	
	#find the f and t indices of each peak
	tind,find = np.where(ispeak)
	
	#create the output array
	n = tind.size
	dtype = [('Power','float32'),('tind','int32'),('find','int32'),
			('t','float32'),('f','float32')]
	out = np.recarray(n,dtype=dtype)
	
	out.Power = P[tind,find+1]
	out.tind = tind
	out.find = find + 1
	out.t = t[tind]
	out.f = f[find + 1]
		
		
	return out
