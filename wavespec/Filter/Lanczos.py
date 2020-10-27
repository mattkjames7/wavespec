import numpy as np


def Lanczos(t0,t,fc,ftype='high'):
	'''
	Calculate the Lanczos-squared function centred upon t0 for an
	arbitrary set of points in time. 
	
	'''
	
	#find the points within 1.5*(1.0/fc) of t0 to apply the function
	#the rest will be set to zero
	n = t.size
	filt = np.zeros(n,dtype='float32')
	pc = 1.0/fc
	use = np.where((t > t0-1.5*pc) & (t < t0+1.5*pc))[0]
	
	
	#get dt
	dt = (t[use] - t0)/(1.5*pc)
	
	#calculate the function
	pidt2 = np.pi*dt
	filt[use] = (np.sin(pidt2)/pidt2)**2
	filt[use] = filt[use]*np.sin(pidt2*3)/(pidt2*3)
	
	#check if any of the dt's are 0
	dt0 = np.where(pidt2 == 0)[0]
	if dt0.size > 0:
		filt[use[dt0]] = 1.0

	#calculate the normalization factor (integrate)
	#tr = t[use[-1]] - t[use[0]]
	norm = np.trapz(filt[use])
	#norm = np.trapz(filt[use],(use.size-1)*(t[use]-t[use[0]])/tr)
	if ftype == 'low':
		filt = filt/norm
	else:
		ic = np.where(t == t0)[0]
		filt = -filt/norm
		filt[ic] = 1.0



	return filt
	
