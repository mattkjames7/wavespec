import numpy as np
import matplotlib.pyplot as plt
from ..Spectrogram.PlotSpectrogram import PlotSpectrogram

def TestSpectrogram():
	
	#define the time axis
	t = np.arange(1000.0)*1.0
	t_diff = np.random.randn(1000)*0.25
	tls = t + t_diff
	tls.sort()
	
	
	#frequency (variable)
	f0 = 0.4
	f1 = 0.03
	ff = 0.005
	f = (f1-f0)*(0.5*np.cos(2*np.pi*ff*t) + 1.0) + f0
	fls = (f1-f0)*(0.5*np.cos(2*np.pi*ff*tls) + 1.0) + f0
	 
	
	#now create the time series
	y = np.cos(2*np.pi*f0*t)*np.cos(2*np.pi*f1*t)
	yls = np.cos(2*np.pi*f0*tls)*np.cos(2*np.pi*f1*tls)
	
	#add a couple of bad data gaps
	y[250:300] = np.nan
	y[550:560] = np.nan
	
	yls[250:300] = np.nan
	yls[550:560] = np.nan
	
	fig = plt
	fig.figure(figsize=(8,11))
	ax0 = plt.subplot2grid((4,1),(0,0))
	ax0.plot(t,y,color=[1.0,0.0,0.0])
	ax0.plot(tls,yls,color=[1.0,0.5,0.0])
	ax0.set_xlim(0,1000)
	fig.subplots_adjust(right=0.9,hspace=0.4,top=0.95,bottom=0.05)
	ax1,_,LenW,freq,_ = PlotSpectrogram(t,y,90.0,6.0,Method='FFT',fig=fig,maps=[1,4,0,1],OneSided=False)
	ax1.set_title('FFT Spectrogram')
	ax1.set_xlim(0,1000)
	ax2,Nw,LenW,Freq,Spec = PlotSpectrogram(t,y,90.0,6.0,Freq=freq[:-1],Method='LS',fig=fig,maps=[1,4,0,2],OneSided=False)
	ax2.set_title('LS Spectrogram')
	ax2.set_xlim(0,1000)
	
	ax3,Nw,LenW,Freq,Spec = PlotSpectrogram(tls,yls,90.0,6.0,Freq=freq[:-1],Method='LS',fig=fig,maps=[1,4,0,3],OneSided=False)
	ax3.set_title('LS Spectrogram (Irregular Data)')
	ax3.set_xlim(0,1000)
	
