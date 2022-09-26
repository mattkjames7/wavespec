import numpy as np
import matplotlib.pyplot as plt
from .Spectrogram import Spectrogram
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ..Tools.DetectGaps import DetectGaps
from ..Tools.UTPlotLabel import UTPlotLabel
import DateTimeTools as TT
from ..Tools.mode import mode
from ..Spectrogram.SpectrogramPlotter import SpectrogramPlotter

#t,v,wind,slip
def PlotSpectrogram(*args,**kwargs):
	'''
	Plots a spectrogram by calling the "Spectrogram" routine which
	creates a spectogram using a sliding window.
	
	Inputs
	======
	t : float
		time array in seconds
	v : float
		array of values the same length as t (when using crossphase
		method, set to a tuple or list containing 2 time series arrays)
	wind : float
		sliding window length in seconds
	slip : float
		difference in time between the start of one window and the 
		next - when slip < wind, each window will have an overlap,
		when slip > wind, there will be gaps where some data will be 
		unused and when slip == wind, each window is adjacent in time.
	Freq : float
		a list of frequencies (Hz) to solve for - only does anything
		when using L-S
	Method : str
		'FFT'|'LS'|'CP-FFT'|'CP-LS' - where CP = crossphase,
		FFT = fast Fourier transform, LS = Lomb-Scargle.
	WindowFunction : str
		Select a window function to apply to the data before 
		the transform, the options are: 'none','cosine-bell','hamming',
		'triangle','welch','blackman','nuttall','blackman-nuttall',
		'flat-top','cosine','gaussian'			
	Param : float
		This parameter is used to alter some of the window functions
		(see WindowFunctions.py).
	Detrend : bool|int
		This will linearly detrend each time window before the 
		transform.
	FindGaps : bool
		This tells the routine to scan for data gaps, when set
		to False - the data are assumed to be perfect.
	GoodData : bool
		This can be set to a boolean array which tells the DetectGaps
		function which data points are good (True) or bad (False),
		if set to None, then any non-finite data is assumed to be bad.
	Quiet : bool
		When set to True, the function produces no stdout output; when
		False, stdout shows the progress.
	fig : object
		This should be set to an instance of matplotlib.pyplot if 
		plotting on an existing figure - if not, then a new figure will
		be created
	maps : list
		This defines the subplot mapping on the figure - 
		[xmaps,ymaps,xmap,ymap], where:
		xmaps - the total number of subplots in x direction
		ymaps - the total number of subplots in y direction
		xmap - the x position of the current subplot, where the starting 
				xmap is 0 which corresponds to the left-most subplot
		ymap - y position of current subplot, 0 being the upper-most 
				plot with increasing integer values moving down the page
	PlotType : str
		'Pow'|'Pha'|'Amp'|'Real'|'Imag'
	scale : float
		Colour scale limits, default is to scale based on the minimum 
		and maximum of the plotted parameter
	zlog : bool
		Set to True for a logarithmic color scale
	TimeAxisUnits : str
		Units of the time axis 's'|'h'|'hh:mm'|'hh:mm:ss'
	FreqAxisUnits : str
		'Hz'|'mHz' - units to use along frequency axis
	Threshold : float
		If set to a value above 0, then all values which 
		correspond to frequencies where the amplitude is less than
		Threshold are set to 0, effectively removing noise from the
		spectra.
	Fudge : bool	(LS Only!)
		This applies a fudge for when f == Nyquist frequency, because
		small floating point numbers have relatively large errors.
		This should only be needed if intending to reproduce a
		two-sided FFT (also, if doing this then divide A by 2 and P 
		by 4).
	OneSided : bool(FFT Only!)
		This should be set to remove the negative frequencies in
		the second half of the spectra. In doing so, the amplitudes
		are doubled and the powers are quadrupled.	
	
	Returns
	=======
	ax : object
		pyploy.Axes instance
	Freq : float
		Array of frequencies in Hz.
	Spec : numpy.recarray 
		Stores the output of the transform under the following fields:
			Tspec : Time in seconds of the middle of each time window
			Pow : Power at each frequency in each window, shape (Nw,LenW)
			Pha : Phase at each frequency in each window, shape (Nw,LenW)
			Amp : Amplitude at each frequency in each window, shape (Nw,LenW)
			Real : Real component at each frequency in each window, shape (Nw,LenW)
			Imag : Imaginary component at each frequency in each window, shape (Nw,LenW)
	'''	

	
	fig = kwargs.get('fig',None)
	maps = kwargs.get('maps',[1,1,0,0])
	PlotType = kwargs.get('PlotType','Pow')
	scale = kwargs.get('scale',None)
	zlog = kwargs.get('zlog',False)
	TimeAxisUnits = kwargs.get('TimeAxisUnits','s')
	FreqAxisUnits = kwargs.get('FreqAxisUnits','Hz')
	nox = kwargs.get('nox',False)
	cmap = kwargs.get('cmap','gnuplot')
	
	if len(args) == 2:
		Freq,Spec = args
		Nw = Spec.size
		Nf = Freq.size - 1
	else:
		t,v,wind,slip = args
		Nw,Freq,Spec = Spectrogram(t,v,wind,slip,**kwargs)
		Nf = Freq.size - 1

	#select the parameter to plot
	#if not PlotType in ['Pow','Pha','Amp','Real','Imag']: 	
		#print('PlotType "{:s}" not recognised - defaulting to "Pow"'.format(PlotType))
		#PlotType = 'Pow'
	if PlotType == 'Real':
		S = Spec.Comp.real()
	elif PlotType == 'Imag':
		S = Spec.Comp.imag()
	else:
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
	dt = mode(ts[1:] - ts[:-1])/2.0
	
	
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

	ax = SpectrogramPlotter(ngd,T0,T1,ts,f,S,fig=fig,maps=maps,zlog=zlog,
									scale=scale,cmap=cmap,zlabel=zlabel)
	
		
	#axis labels
	ax.set_xlabel(xlabel)
	ax.set_ylabel('$f$ ('+FreqAxisUnits+')')
		
	#sort the time axis out
	if nox:
		ax.xaxis.set_visible(False)
	else:
		if TimeAxisUnits in ['hh:mm','hh:mm:ss']:
			TT.DTPlotLabel(ax)
			ax.set_xlabel('UT')	

			
	return ax,Freq,Spec
