import numpy as np
import matplotlib.pyplot as plt

def RemoveStep(t,xin,step,Order=2,nInd=5,Gap=1):
	'''
	Remove any annoying steps in data (this requires an integer array 
	which says where the steps are)
	
	'''
	
	#firstly we need to work out where the steps are
	x = np.copy(xin)
	si0 = np.where(step == 1)[0]
	si1 = np.where(step == 2)[0]

	if si0.size == 0 or si1.size == 0:
		return x

	#now work out where the continuous bits of good data are
	if si0[0] < si1[0]:
		g0 = np.append(0,si1)
	else:
		g0 = si1
	
	if si1[-1] > si0[-1]:
		g1 = np.append(si0,x.size-1)
	else:
		g1 = si0
		
	#plt.figure()
	#plt.plot(t,x,color='blue')
	ng = g0.size
	#loop through each pair of good bits 
	for i in range(0,ng-1):
		#get last nInd points of the first good bit and the first nInd
		#points of the next good bit...
		#first bit
		a0 = np.max([g0[i],g1[i] - nInd - Gap])
		a1 = g1[i] + 1 - Gap
		ta = t[a0:a1]
		xa = x[a0:a1]
		#second bit
		b0 = g0[i+1] + Gap
		b1 = np.min([g1[i+1]+1+Gap,g0[i+1]+1+nInd+Gap])
		tb = t[b0:b1]
		xb = x[b0:b1]
		
		#plt.scatter(ta,xa,color='red')
		#plt.scatter(tb,xb,color='orange')
		
		
		#get the polyfits
		Oa = np.min([Order,xa.size-1])
		Ob = np.min([Order,xb.size-1])
		if Oa == 0:
			x0 = xa[-1]
		else:
			pfa = np.polyfit(ta,xa,Oa)
			Pa = np.poly1d(pfa)
			x0 = Pa(tb[0])
			#plt.plot(ta,Pa(ta),color='red')
		if Ob == 0:
			x1 = xb[0]
		else:
			pfb = np.polyfit(tb,xb,Ob)
			Pb = np.poly1d(pfb)
			x1 = Pb(tb[0])
			#plt.plot(tb,Pb(tb),color='orange')
			
		#subtract this from all of the next points
		dx = x1 - x0
		x[b0:] -= dx

		#plt.plot(t[b0:],x[b0:],color='green')

		#assess whether each point from a1 - Gap to b0 + Gap [a1-Gap:b0]
		#is closer to original or new points
		tg = t[a1:b0]
		xg = x[a1:b0]
		#plt.scatter(tg,xg,color='pink')
		if tg.size > 0:
			for j in range(0,tg.size):
				if Oa == 0:
					da = np.abs(xg[j] - x0)
				else:
					da = np.abs(xg[j] - Pa(tg[j]))
				if Ob == 0:
					db = np.abs(xg[j] - x1)
				else:
					db = np.abs(xg[j] - Pb(tg[j]))
				if da > db:
					#this means that we should apply -dx
					x[a1 + j] -= dx
	#plt.plot(t,x,color='black')
	return x
		
