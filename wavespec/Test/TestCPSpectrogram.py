import numpy as np
import matplotlib.pyplot as plt
from ..Spectrogram.PlotSpectrogram import PlotSpectrogram

def TestCPSpectrogram():
	
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
	 
	dphi = 5.0*np.pi/180.0
	#now create the time series
	y0 = np.cos(2*np.pi*f0*t)*np.cos(2*np.pi*f1*t)
	yls0 = np.cos(2*np.pi*f0*tls)*np.cos(2*np.pi*f1*tls)
	y1 = np.cos(2*np.pi*f0*t)*np.cos(2*np.pi*f1*t+dphi)
	yls1 = np.cos(2*np.pi*f0*tls)*np.cos(2*np.pi*f1*tls+dphi)
	
	#add a couple of bad data gaps
	y0[250:300] = np.nan
	y0[550:560] = np.nan
	
	y1[250:300] = np.nan
	y1[550:560] = np.nan
	
	yls0[250:300] = np.nan
	yls0[550:560] = np.nan
	
	yls1[250:300] = np.nan
	yls1[550:560] = np.nan
	
	y = (y0,y1)
	yls = (yls0,yls1)
	
	fig = plt
	fig.figure(figsize=(8,8))
	ax0 = plt.subplot2grid((3,1),(0,0))
	ax0.plot(t,y[0],color=[1.0,0.0,0.0])
	ax0.plot(t,y[1],color=[1.0,0.0,0.0])
	ax0.plot(tls,yls[0],color=[1.0,0.5,0.0])
	ax0.plot(tls,yls[1],color=[1.0,0.5,0.0])
	ax0.set_xlim(0,1000)
	fig.subplots_adjust(right=0.8,hspace=0.3)
	ax1,_,_,freq,_ = PlotSpectrogram(t,y,90.0,6.0,Method='CP-FFT',fig=fig,maps=[1,3,0,1])
	ax1.set_title('FFT Spectrogram')
	ax1.set_xlim(0,1000)
	ax2,Nw,LenW,Freq,Spec = PlotSpectrogram(tls,yls,90.0,6.0,Freq=freq,Method='CP-LS',fig=fig,maps=[1,3,0,2])
	ax2.set_title('LS Spectrogram')
	ax2.set_xlim(0,1000)
	
