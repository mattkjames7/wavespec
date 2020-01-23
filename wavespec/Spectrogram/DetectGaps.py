import numpy as np


def DetectGaps(v,GoodData=None):
	'''
	Detects gaps in a time series either by looking for inf/nan in the
	input array, or using a supplimentary boolean array.
	
	Inputs
	======
	v:	Input time series.
	GoodData: Optional - When set to None, it is ignored, gaps in v are
		assumed to be not finite, otherwise set it to a Boolean array
		the same length as v, where True is good and False is bad.
		
	Returns
	=======
	ngd:	The number of good sections of data.
	UTi0:	Start index of a good bit of data.
	UTi1:	End index of a good section of data.
	
	'''
	Tlen = v.size
	if GoodData is None:
		good = np.zeros(Tlen,dtype='bool')
		gd = np.where(np.isfinite(v))[0]
		good[gd] = True
	else:
		good = GoodData
	st = -1
	ngd = 0
	for i in range(0,Tlen):
		if good[i]:
			if st == -1:
				st = i
		else:
			if st != -1:
				if ngd == 0:
					UTi0 = np.array([st])
					UTi1 = np.array([i-1])
				else:
					UTi0 = np.append(UTi0,st)
					UTi1 = np.append(UTi1,i-1)	
				st = -1
				ngd += 1
	if st != -1:
		if ngd == 0:
			UTi0 = np.array([st])
			UTi1 = np.array([i-1])
		else:
			UTi0 = np.append(UTi0,st)
			UTi1 = np.append(UTi1,i)	
		st =-1
		ngd += 1
	
		
	if ngd == 0:
		UTi0 = np.array([])
		UTi1 = np.array([])
		
	
	return ngd,UTi0,UTi1
