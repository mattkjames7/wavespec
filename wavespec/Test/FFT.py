'''
Test the FFT functions
'''

import numpy as np
from .. import Fourier
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ..Tools.mode import mode
from ..Tools.DetectGaps import DetectGaps
from ..Spectrogram.SpectrogramPlotter import SpectrogramPlotter

def Spectrum():
	
	#pick two frequencies
	f0 = 0.002
	f1 = 0.005
	
	#amplitudes
	A0 = 2.0
	A1 = 1.5

	#phases
	p0 = np.pi/2.0
	p1 = 0.0
	
	#time series
	t = np.arange(1000.0)
	v0 = A0*np.cos(2*np.pi*f0*t + p0)
	v1 = A1*np.cos(2*np.pi*f1*t + p1)
	v = v0 + v1
	
	#spectrum
	power,A,phase,fr,fi,freq = Fourier.FFT(t,v,OneSided=True)
	
	#figure
	fig = plt
	fig.figure(figsize=(8,11))
	
	ax0 = fig.subplot2grid((2,1),(0,0))
	ax1 = fig.subplot2grid((2,1),(1,0))
	
	ax0.plot(t,v0,color='red',linestyle='--')
	ax0.plot(t,v1,color='orange',linestyle='--')
	ax0.plot(t,v,color='black',linestyle='-')
	
	ax1.plot(freq[:-1]*1000.0,power,color='blue')
	
	ax0.set_xlabel('$t$ (s)')
	
	ax1.set_ylabel('Power')
	ax1.set_xlabel('Frequency (mHz)')
	
	fmx = np.min([freq.max(),1.5*np.max([f0,f1])])
	ax1.set_xlim(0,fmx*1000)
	
	
def Spectrogram():
	
	#pick two frequencies
	f0 = 0.002
	f1 = 0.005
	
	#amplitudes
	A0 = 2.0
	A1 = 1.5

	#phases
	p0 = np.pi/2.0
	p1 = 0.0
	
	#time series
	t = np.arange(10800.0)
	v0 = A0*np.cos(2*np.pi*f0*t + p0)
	v1 = A1*np.cos(2*np.pi*f1*t + p1)
	v = v0 + v1
	
	wind = 1800
	slip = 200

	#figure
	fig = plt
	fig.figure(figsize=(8,11))
	
	ax0 = fig.subplot2grid((2,1),(0,0))
	ax1 = fig.subplot2grid((2,1),(1,0))

	ax0.plot(t,v0,color='red',linestyle='--')
	ax0.plot(t,v1,color='orange',linestyle='--')
	ax0.plot(t,v,color='black',linestyle='-')
	ax0.set_xlabel('Time (s)')
		
	#spectrogram
	ax1,Freq,Spec = Fourier.PlotSpectrogram(t,v,wind,slip,FreqAxisUnits='mHz',fig=fig,maps=[1,2,0,1],TimeAxisUnits='s')
	fmx = np.min([Freq.max(),1.5*np.max([f0,f1])])
	ax1.set_ylim(0,fmx*1000)
	
def Spectrogram2():
	
	#pick two frequencies
	f0 = 0.002
	f1 = 0.005
	
	#amplitudes
	A0 = 2.0
	A1 = 1.5

	#phases
	p0 = np.pi/2.0
	p1 = 0.0
	
	#time series
	t = np.arange(10800.0)
	v0 = A0*np.cos(2*np.pi*f0*t + p0)
	v1 = A1*np.cos(2*np.pi*f1*t + p1)
	v = v0 + v1
	
	wind = 1800
	slip = 200

	#figure
	fig = plt
	fig.figure(figsize=(8,11))
	
	ax0 = fig.subplot2grid((2,1),(0,0))
	ax1 = fig.subplot2grid((2,1),(1,0))

	ax0.plot(t,v0,color='red',linestyle='--')
	ax0.plot(t,v1,color='orange',linestyle='--')
	ax0.plot(t,v,color='black',linestyle='-')
	ax0.set_xlabel('Time (s)')
		
	Nw,Freq,Spec = Fourier.Spectrogram(t,v,wind,slip)
		
	#spectrogram
	ax1,Freq,Spec = Fourier.PlotSpectrogram(Freq,Spec,FreqAxisUnits='mHz',fig=fig,maps=[1,2,0,1])
	fmx = np.min([Freq.max(),1.5*np.max([f0,f1])])
	ax1.set_ylim(0,fmx*1000)
	
def Spectrogram3():
	
	#pick two frequencies
	f0 = 0.002
	f1 = 0.005
	
	#amplitudes
	A0 = 2.0
	A1 = 1.5

	#phases
	p0 = np.pi/2.0
	p1 = 0.0
	
	#time series
	t = np.arange(10800.0)
	v0 = A0*np.cos(2*np.pi*f0*t + p0)
	v1 = A1*np.cos(2*np.pi*f1*t + p1)
	v = v0 + v1
	
	wind = 1800
	slip = 200

	#figure
	fig = plt
	fig.figure(figsize=(8,11))
	
	ax0 = fig.subplot2grid((2,1),(0,0))
	ax1 = fig.subplot2grid((2,1),(1,0))

	ax0.plot(t,v0,color='red',linestyle='--')
	ax0.plot(t,v1,color='orange',linestyle='--')
	ax0.plot(t,v,color='black',linestyle='-')
	ax0.set_xlabel('Time (s)')
		
	Nw,Freq,Spec = Fourier.Spectrogram(t,v,wind,slip)
		
	#spectrogram
	ax1 = SpectrogramPlotter(Spec.Tspec,Freq*1000,Spec.Pow,fig=fig,maps=[1,2,0,1])
	fmx = np.min([Freq.max(),1.5*np.max([f0,f1])])
	ax1.set_ylim(0,fmx*1000)
	

def Spectrogram3D():


	#pick some frequencies
	fx0 = 0.002
	fx1 = 0.007
	fy0 = 0.007
	fy1 = 0.010
	
	#amplitudes
	A0 = 2.0
	A1 = 1.5

	#phases
	p0 = np.pi/2.0
	p1 = 0.0
	
	#time series
	t = np.arange(10800.0)
	x0 = A0*np.cos(2*np.pi*fx0*t + p0)
	x1 = A1*np.cos(2*np.pi*fx1*t + p1)
	x = x0 + x1
	y0 = A0*np.cos(2*np.pi*fy0*t + p0)
	y1 = A1*np.cos(2*np.pi*fy1*t + p1)
	y = y0 + y1
	z = np.zeros(t.size,dtype='float32')
	
	#spectrogram
	wind = 1800
	slip = 200
	Nw,Freq,Spec = Fourier.Spectrogram3D(t,x,y,z,wind,slip,CombineComps=True)
	Nf = Freq.size - 1
	S = Spec.xyPow
	f = Freq[:Nf+1]*1000.0
	ts = Spec.Tspec
	xlabel = 'Time (s)'
	dt = mode(ts[1:] - ts[:-1])/2.0

	scale = [np.nanmin(S),np.nanmax(S)]
	norm = colors.Normalize(vmin=scale[0],vmax=scale[1])	
	cmap = plt.cm.get_cmap('gnuplot')

	#find gaps
	gaps = np.where(np.isfinite(S[:,1]) == False)[0]
	ngd,T0,T1 = DetectGaps(S[:,1])


	#figure
	fig = plt
	fig.figure(figsize=(8,11))
	
	ax0 = fig.subplot2grid((2,1),(0,0))
	ax1 = fig.subplot2grid((2,1),(1,0))
	ax0.plot(t,x,color='red')
	ax0.plot(t,y,color='orange')
	
	
	sm = None
	for i in range(0,ngd):
		#select the good portion of the 
		use = np.arange(T0[i],T1[i]+1)
		tax = np.append(ts[use]-dt,ts[use[-1]]+dt)
		Stmp = S[use]
		
		
		#mesh the axes
		tm,fm = np.meshgrid(tax,f)
		#plot the section
		sm = ax1.pcolormesh(tm.T,fm.T,Stmp,cmap=cmap,norm=norm)


	#colour bar
	fig.subplots_adjust(right=0.8)
	box = ax1.get_position()
	if not sm is None:
		cax = plt.axes([0.05*box.width + box.x1,box.y0+0.1*box.height,box.width*0.025,box.height*0.8])
		cbar = fig.colorbar(sm,cax=cax)
		cbar.set_label('Power')
		
	#axis labels
	ax1.set_xlabel(xlabel)
	ax1.set_ylabel('$f$ (mHz)')


	fmx = np.min([Freq.max(),1.5*np.max([fx0,fx1,fy0,fy1])])
	ax1.set_ylim(0,fmx*1000)

	return Spec

def Spectrogram3D2():


	#pick some frequencies
	fx0 = 0.007
	fx1 = 0.002
	fy0 = 0.007
	fy1 = 0.010
	
	#amplitudes
	A0 = 2.0
	A1 = 1.5

	#phases
	px0 = np.append(np.linspace(0.0,2*np.pi,5400),np.zeros(5400))
	py0 = np.append(np.linspace(0.0,2*np.pi,5400)-np.pi/2.0,np.linspace(0.0,2*np.pi,5400))
	p1 = 0.0
	
	#time series
	t = np.arange(10800.0)
	x0 = A0*np.cos(2*np.pi*fx0*t + px0)
	x1 = A1*np.cos(2*np.pi*fx1*t + p1)
	x = x0 + x1
	y0 = A0*np.cos(2*np.pi*fy0*t + py0)
	y1 = A1*np.cos(2*np.pi*fy1*t + p1)
	y = y0 + y1
	z = np.zeros(t.size,dtype='float32')
	
	#spectrogram
	wind = 1800
	slip = 200
	Nw,Freq,Spec = Fourier.Spectrogram3D(t,x,y,z,wind,slip,CombineComps=True)
	# Nf = Freq.size - 1
	# S = Spec.xyPow
	# f = Freq[:Nf+1]*1000.0
	# ts = Spec.Tspec
	xlabel = 'Time (s)'
	# dt = mode(ts[1:] - ts[:-1])/2.0

	# scale = [np.nanmin(S),np.nanmax(S)]
	# norm = colors.Normalize(vmin=scale[0],vmax=scale[1])	
	# cmap = plt.cm.get_cmap('gnuplot')

	# #find gaps
	# gaps = np.where(np.isfinite(S[:,1]) == False)[0]
	# ngd,T0,T1 = DetectGaps(S[:,1])


	#figure
	fig = plt
	fig.figure(figsize=(11,8))
	plt.subplots_adjust(left=0.05,right=1)
	
	ax0 = fig.subplot2grid((2,2),(0,0))
	ax1 = fig.subplot2grid((2,2),(1,0))
	ax2 = fig.subplot2grid((2,2),(0,1))
	ax3 = fig.subplot2grid((2,2),(1,1))
	ax0.plot(t,x,color='red')
	ax0.plot(t,y,color='orange')
	
	ax1 = SpectrogramPlotter(Spec.Tspec,Freq*1000,Spec.xyPow,fig=fig,maps=[2,2,0,1],zlabel='xyPower')
	ax2 = SpectrogramPlotter(Spec.Tspec,Freq*1000,Spec.xPow,fig=fig,maps=[2,2,1,0],zlabel='xPower')
	ax3 = SpectrogramPlotter(Spec.Tspec,Freq*1000,Spec.yPow,fig=fig,maps=[2,2,1,1],zlabel='yPower')
	
	# sm = None
	# for i in range(0,ngd):
		# #select the good portion of the 
		# use = np.arange(T0[i],T1[i]+1)
		# tax = np.append(ts[use]-dt,ts[use[-1]]+dt)
		# Stmp = S[use]
		
		
		# #mesh the axes
		# tm,fm = np.meshgrid(tax,f)
		# #plot the section
		# sm = ax1.pcolormesh(tm.T,fm.T,Stmp,cmap=cmap,norm=norm)


	# #colour bar
	# fig.subplots_adjust(right=0.8)
	# box = ax1.get_position()
	# if not sm is None:
		# cax = plt.axes([0.05*box.width + box.x1,box.y0+0.1*box.height,box.width*0.025,box.height*0.8])
		# cbar = fig.colorbar(sm,cax=cax)
		# cbar.set_label('Power')
		
	#axis labels
	ax1.set_xlabel(xlabel)
	ax1.set_ylabel('$f$ (mHz)')
	ax2.set_xlabel(xlabel)
	ax2.set_ylabel('$f$ (mHz)')
	ax3.set_xlabel(xlabel)
	ax3.set_ylabel('$f$ (mHz)')


	fmx = np.min([Freq.max(),1.5*np.max([fx0,fx1,fy0,fy1])])
	ax1.set_ylim(0,fmx*1000)
	ax2.set_ylim(0,fmx*1000)
	ax3.set_ylim(0,fmx*1000)

	#plt.tight_layout()
	#return Spec
