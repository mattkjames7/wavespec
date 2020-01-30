import numpy as np
import matplotlib.pyplot as plt
from .UTPlotLabel import UTPlotLabel

def PlotPolarization(t,Ax,Ay,Px,Py,Dir,fig=None,maps=[1,1,0,0],Multiplier=1.0,nox=False,trange=None,TimeAxisUnits='s'):
	
	if fig == None:
		fig = plt
		fig.figure()
		
	ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
	bbsize = ax.bbox.size

	#scale the time axis
	if TimeAxisUnits == 'h':
		ts = t/3600.0
		xlabel = 'Time (h)'
	elif TimeAxisUnits in ['hh:mm','hh:mm:ss']:
		ts = t/3600.0
		xlabel = 'Time'
	else:
		ts = t
		xlabel = 'Time (s)'

	if trange is None:
		trange = [ts[0],ts[-1]]
	tr = trange[1] - trange[0]
	yr = 0.5*tr*bbsize[1]/bbsize[0]
	
	ax.axis([trange[0],trange[1],-yr,yr])
	
	Mx = np.nanmax([Ax,Ay])
	MaxAmp = yr*Multiplier
	axx = MaxAmp*Ax/Mx
	ayy = MaxAmp*Ay/Mx

	n = np.size(t)
	a = np.arange(360.0)*np.pi/180.0
	for i in range(0,n):
		x = ts[i] + axx[i]*np.cos(a + Px[i])
		y = ayy[i]*np.cos(a + Py[i])
		if Dir[i] > 0:
			col = [1.0,0.0,0.0]
		else:
			col = [0.0,1.0,0.0]
		
		fig.plot(x,y,color=col,linewidth=1.0)
	
	R = ax.axis()

	fig.annotate('RH',xy=(0.99*(R[1]-R[0])+R[0],R[2]+0.9*(R[3]-R[2])),color=[1.0,0.0,0.0],ha='right')	
	fig.annotate('LH',xy=(0.99*(R[1]-R[0])+R[0],R[2]+0.8*(R[3]-R[2])),color=[0.0,1.0,0.0],ha='right')
	if nox:
		ax.xaxis.set_visible(False)
	else:
		if TimeAxisUnits in ['hh:mm','hh:mm:ss']:
			UTPlotLabel(ax,axis='x',seconds=(TimeAxisUnits == 'hh:mm:ss'))
			ax.set_xlabel(xlabel)
	lbl = ax.get_yticklabels()
	ax.set_yticklabels(['']*np.size(lbl))
	
	return ax
		
