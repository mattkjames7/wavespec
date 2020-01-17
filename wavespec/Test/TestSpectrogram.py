import numpy as np
import matplotlib.pyplot as plt
from ..Spectrogram.PlotSpectrogram import PlotSpectrogram

def TestSpectrogram():
	
	#define the time axis
	t = np.arange(1000.0)*1.0
	
	#frequency (variable)
	f0 = 0.4
	f1 = 0.03
	ff = 0.005
	f = (f1-f0)*(0.5*np.cos(2*np.pi*ff*t) + 1.0) + f0

	
	#now create the time series
	y = np.cos(2*np.pi*f0*t)*np.cos(2*np.pi*f1*t)
	
	#add a couple of bad data gaps
	y[250:300] = np.nan
	y[550:560] = np.nan
	
	fig = plt
	fig.figure(figsize=(8,8))
	ax0 = plt.subplot2grid((3,1),(0,0))
	ax0.plot(t,y)
	ax0.set_xlim(0,1000)
	fig.subplots_adjust(right=0.8,hspace=0.3)
	ax1,_,_,freq,_ = PlotSpectrogram(t,y,90.0,6.0,Method='FFT',fig=fig,maps=[1,3,0,1])
	ax1.set_title('FFT Spectrogram')
	ax1.set_xlim(0,1000)
	ax2,Nw,LenW,Freq,Spec = PlotSpectrogram(t,y,90.0,6.0,Freq=freq,Method='LS',fig=fig,maps=[1,3,0,2])
	ax2.set_title('LS Spectrogram')
	ax2.set_xlim(0,1000)
	
#	fig.tight_layout()
	return Nw,LenW,Freq,Spec
