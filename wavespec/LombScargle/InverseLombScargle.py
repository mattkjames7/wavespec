import numpy as np


def InverseLombScargle(t,f,A,phi):
	'''
	A highly experimental way of reversing the L-S periodogram to 
	reproduce the data. It does not work.
	
	
	Inputs
	======
	t : float
		Times at which to reproduce the time series
	f : float
		Frequencies at which the periodogram was created
	A : float
		Amplitudes from the periodogram
	phi : float
		Phases produces by the periodogram
	NOTE: Reversing a periodogram which had a window function applied to
	will have odd effects.
	
	'''
	
	#create an output array
	x = np.zeros(t.size)
	#tm = 0.5*(t[1] + t[0])
	
	#now do a sum of cosines
	w = 2*np.pi*f
	for i in range(0,w.size):
		if w[i] == 0:
			x[:] += A[i]
		else:
			x += A[i]*np.cos(w[i]*(t-t[0]) - phi[i])
	
	return x
