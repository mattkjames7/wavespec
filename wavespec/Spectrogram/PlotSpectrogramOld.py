import numpy as np
import matplotlib.pyplot as plt
from .Spectrogram import Spectrogram
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from .DetectGaps import DetectGaps
from ..Tools.UTPlotLabel import UTPlotLabel
import DateTimeTools as TT

from .. import Fourier
from .. import LombScargle
from .. import CrossPhase

def _mode(x):
	u,c = np.unique(x,return_counts=True)
	return u[c.argmax()]
	


def PlotSpectrogram(t,v,wind,slip,Freq=None,Method='FFT',WindowFunction=None,
					Param=None,Detrend=True,FindGaps=True,GoodData=None,
					Quiet=True,LenW=None,fig=None,maps=[1,1,0,0],PlotType='Pow',
					scale=None,zlog=False,TimeAxisUnits='s',FreqAxisUnits='Hz',
					nox=False,Threshold=0.0,Fudge=False,OneSided=True):
	'''
	Plots a spectrogram by calling the "Spectrogram" routine which
	creates a spectogram using a sliding window.
	
	Inputs
	======
		t : time array in seconds
		v : array of values the same length as t (when using crossphase
			method, set to a tuple or list containing 2 time series arrays)
	 wind : sliding window length in seconds
	 slip : difference in time between the start of one window and the 
			next - when slip < wind, each window will have an overlap,
			when slip > wind, there will be gaps where some data will be 
			unused and when slip == wind, each window is adjacent in time.
	 Freq : a list of frequencies (Hz) to solve for - only does anything
			when using L-S
	Method : 'FFT'|'LS'|'CP-FFT'|'CP-LS' - where CP = crossphase,
			FFT = fast Fourier transform, LS = Lomb-Scargle.
	WindowFunction : Select a window function to apply to the data before 
			the transform, the options are: 'none','cosine-bell','hamming',
			'triangle','welch','blackman','nuttall','blackman-nuttall',
			'flat-top','cosine','gaussian'			
	Param : This parameter is used to alter some of the window functions
			(see WindowFunctions.py).
	Detrend : This will linearly detrend each time window before the 
			transform.
	FindGaps : This tells the routine to scan for data gaps, when set
			to False - the data are assumed to be perfect.
	GoodData : This can be set to a boolean array which tells the DetectGaps
			function which data points are good (True) or bad (False),
			if set to None, then any non-finite data is assumed to be bad.
	Quiet : When set to True, the function produces no stdout output; when
			False, stdout shows the progress.
	LenW : This can be set to an integer value in order to for a specific
			window length (the number of elements, as opposed to the length
			in time defined using wind)
	fig : This should be set to an instance of matplotlib.pyplot if 
			plotting on an existing figure - if not, then a new figure will
			be created
	maps : This defines the subplot mapping on the figure - 
			[xmaps,ymaps,xmap,ymap], where:
			xmaps - the total number of subplots in x direction
			ymaps - the total number of subplots in y direction
			xmap - the x position of the current subplot, where the starting 
					xmap is 0 which corresponds to the left-most subplot
			ymap - y position of current subplot, 0 being the upper-most 
					plot with increasing integer values moving down the page
	PlotType : 'Pow'|'Pha'|'Amp'|'Real'|'Imag'
	scale : Colour scale limits, default is to scale based on the minimum 
			and maximum of the plotted parameter
	zlog : Set to True for a logarithmic color scale
	TimeAxisUnits : Units of the time axis 's'|'h'|'hh:mm'|'hh:mm:ss'
	FreqAxisUnits : 'Hz'|'mHz' - units to use along frequency axis
	Threshold:	If set to a value above 0, then all values which 
			correspond to frequencies where the amplitude is less than
			Threshold are set to 0, effectively removing noise from the
			spectra.
	Fudge:	(LS Only!)
			This applies a fudge for when f == Nyquist frequency, because
			small floating point numbers have relatively large errors.
			This should only be needed if intending to reproduce a
			two-sided FFT (also, if doing this then divide A by 2 and P 
			by 4).
	OneSided: (FFT Only!)
			This should be set to remove the negative frequencies in
			the second half of the spectra. In doing so, the amplitudes
			are doubled and the powers are quadrupled.	
	
	Returns
	=======
	Nw : Total number of time windows in the output array
	LenW : Length of a time window (number of elements)
	Freq : Array of frequencies in Hz.
	numpy.recarray : 
			Stores the output of the transform under the following fields:
				Tspec : Time in seconds of the middle of each time window
				Pow : Power at each frequency in each window, shape (Nw,LenW)
				Pha : Phase at each frequency in each window, shape (Nw,LenW)
				Amp : Amplitude at each frequency in each window, shape (Nw,LenW)
				Real : Real component at each frequency in each window, shape (Nw,LenW)
				Imag : Imaginary component at each frequency in each window, shape (Nw,LenW)
	'''	
	
	
	Nw,LenW,Freq,Spec = Spectrogram(t,v,wind,slip,Freq,Method,WindowFunction,Param,Detrend,FindGaps,GoodData,Quiet,LenW,Threshold=Threshold,Fudge=Fudge,OneSided=OneSided)
	Nf = Freq.size - 1

	#select the parameter to plot
	if not PlotType in ['Pow','Pha','Amp','Real','Imag']: 	
		print('PlotType "{:s}" not recognised - defaulting to "Pow"'.format(PlotType))
		PlotType = 'Pow'
	S = Spec[PlotType]
	
	#scale the time axis
	if TimeAxisUnits == 'h':
		ts = Spec.Tspec/3600.0
		xlabel = 'Time (h)'
	elif TimeAxisUnits in ['hh:mm','hh:mm:ss']:
		ts = Spec.Tspec/3600.0
		xlabel = 'Time'
	else:
		ts = Spec.Tspec
		xlabel = 'Time (s)'
	dt = _mode(ts[1:] - ts[:-1])/2.0
	
	
	#find gaps
	gaps = np.where(np.isfinite(Spec.Pow[:,1]) == False)[0]
	ngd,T0,T1 = DetectGaps(Spec.Pow[:,1])
	

	
	#set the frequency
	if FreqAxisUnits == 'Hz':
		f = Freq[:Nf+1]
	elif FreqAxisUnits == 'mHz':
		f = Freq[:Nf+1]*1000.0
	else:
		print('Frequency axis units {:s} not recognised, defaulting to "Hz"'.format(FreqAxisUnits))
		f = Freq[:Nf+1]
	
	#set the z (colour) scale
	zunits = { 'Pow' : 'Power',
			   'Pha' : 'Phase',
			   'Amp' : 'Amplitude',
			   'Real' : 'Real Component',
			   'Imag' : 'Imaginary Component'}
	zlabel = zunits[PlotType]
	if scale is None:
		scale = [np.nanmin(S),np.nanmax(S)]
	if zlog:
		if scale == 0.0:
			scale[0] = np.nanmin(S[(S > 0) & np.isfinite(S)])
		norm = colors.LogNorm(vmin=scale[0],vmax=scale[1])

	else:
		norm = colors.Normalize(vmin=scale[0],vmax=scale[1])	
	
	#create the plot
	if fig is None:
		fig = plt
		fig.figure()
	ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
	cmap = plt.cm.get_cmap('gnuplot')
	
	
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

	#colour bar
	fig.subplots_adjust(right=0.8)
	box = ax.get_position()
	if not sm is None:
		cax = plt.axes([0.05*box.width + box.x1,box.y0+0.1*box.height,box.width*0.025,box.height*0.8])
		cbar = fig.colorbar(sm,cax=cax)
		cbar.set_label(zlabel)
		
	#axis labels
	ax.set_xlabel(xlabel)
	ax.set_ylabel('$f$ ('+FreqAxisUnits+')')
		
	#sort the time axis out
	if nox:
		ax.xaxis.set_visible(False)
	else:
		if TimeAxisUnits in ['hh:mm','hh:mm:ss']:
			#UTPlotLabel(ax,axis='x',seconds=(TimeAxisUnits == 'hh:mm:ss'))
			TT.DTPlotLabel(ax)
			ax.set_xlabel('UT')	

			
	return ax,Nw,LenW,Freq,Spec
