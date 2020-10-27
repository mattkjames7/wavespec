import numpy as np


def PolyDetrend(t,x,Order=1):
	'''
	Detrends a time series using a polynomial.
	
	Inputs
	======
	t : float
		Time array
	x : float
		Value array to be detrended
	Order : int
		Order of polynomial to use
		
	Returns
	=======
	y : float
		Detrended x
	'''


	#get the polyfit first
	good = np.isfinite(x)
	if not good.all():
		use = np.where(good)[0]
		p = np.polyfit(t[good],x[good],Order)
	else:
		p = np.polyfit(t,x,Order)
	
	#get the polynomial function
	pf = np.poly1d(p)
	
	#now work out the values of the fit at y
	xp = pf(t)
	
	#subtract from original
	return x - xp
