import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ..Tools.mode import mode

def SpectrogramPlotter(*args,**kwargs):

	fig = kwargs.get('fig',None)
	maps = kwargs.get('maps',[1,1,0,0])
	scale = kwargs.get('scale',None)
	zlog = kwargs.get('zlog',False)
	zlabel = kwargs.get('zlabel','')
	cmap = kwargs.get('cmap','gnuplot')
	

	
	#two options:
	if len(args) == 3:
		#1 - simple taxis,faxis,spec
		ts,f,S = args
		if ts.size == S.shape[0]:
			dt = mode(ts[1:] - ts[:-1])/2.0
			tax = np.append(ts-dt,ts[-1]+dt)
		else:
			tax = ts
		tlim = [tax[0],tax[-1]]

		MultiSpec = False
	else:
		#2 - split into multiple spectra
		ngd,T0,T1,ts,f,S = args
		dt = mode(ts[1:] - ts[:-1])/2.0
		MultiSpec = True
		tlim = [ts[0],ts[-1]+dt]

	#check that f is the right length
	if f.size == S.shape[1]:
		df = f[1] - f[0]
		f = np.append(f,f[-1]+df)


	if scale is None:
		if zlog:
			lS = np.log10(S)
			lS[np.isinf(lS)] = np.nan
			muS = np.nanmean(lS)
			sgS = np.nanstd(lS)
			scale = 10**np.array([muS-1.0*sgS,muS+1.0*sgS])
			if (np.isfinite(scale) == False).any():
				scale=[0.1,10.0]
		else:
			scale = [np.nanmin(S),np.nanmax(S)]
	if zlog:
		if scale[0] == 0.0:
			scale[0] = np.nanmin(S[(S > 0) & np.isfinite(S)])
		norm = colors.LogNorm(vmin=scale[0],vmax=scale[1])
	else:
		norm = colors.Normalize(vmin=scale[0],vmax=scale[1])	
	
	#create the plot
	if fig is None:
		fig = plt
		fig.figure()
	ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
	

	if MultiSpec:
		#loop through each good section
		sm = None
		for i in range(0,ngd):
			#select the good portion of the 
			use = np.arange(T0[i],T1[i]+1)
			tax = np.append(ts[use]-dt,ts[use[-1]]+dt)
			Stmp = S[use]
			
			
			#mesh the axes
			tm,fm = np.meshgrid(tax,f)
			#plot the section
			sm = ax.pcolormesh(tm.T,fm.T,Stmp,cmap=cmap,norm=norm)
	else:
		#mesh the axes
		tm,fm = np.meshgrid(tax,f)
		#plot the section
		sm = ax.pcolormesh(tm.T,fm.T,S,cmap=cmap,norm=norm)		

	
	#colour bar
	fig.subplots_adjust(right=0.8)
	box = ax.get_position()
	if not sm is None:
		cax = plt.axes([0.05*box.width + box.x1,box.y0+0.1*box.height,box.width*0.025,box.height*0.8])
		cbar = fig.colorbar(sm,cax=cax)
		cbar.set_label(zlabel)
	
	ax.set_xlim(tlim)
	return ax
