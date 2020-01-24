import numpy as np
import matplotlib.pyplot as plt
from ..Fourier import FFT
from ..CrossPhase.CrossPhase import CrossPhase

def TestCP(f=0.08,phi=[0.0,50.0],Backend='C++'):
	'''
	Test the Lomb-Scargle routine against the FFT.
	 
	'''
	#create the time series
	A = 1.0
	n = np.size(phi)
	t = np.arange(0,51.0)
	y = np.zeros((n,51),dtype='float32')
	for i in range(0,n):
		y[i] = A*np.cos(2*np.pi*f*t + np.pi*phi[i]/180.0)
	
	mx = np.max(np.abs(y))
	
	fig = plt
	fig.figure(figsize=(8,8))
	ax0 = fig.subplot2grid((2,1),(0,0))
	for i in range(0,n):
		ax0.plot(t,y[i],label='$A_{'+'{:d}'.format(i)+'}$='+'{:3.1f}'.format(A)+' $f_{'+'{:d}'.format(i)+'}$='+'{:4.2f}'.format(f)+r' $\phi_{'+'{:d}'.format(i)+'}$='+'{:5.1f}'.format(phi[i])) 
	ax0.plot(t,y[0],label='$y_0$')
	ax0.plot(t,y[1],label='$y_1$')
	ax0.set_xlabel('$t$')
	ax0.set_ylabel('$y$')
	ax0.plot([0.0,50.0],[0.0,0.0],color=[0.0,0.0,0.0],linestyle=':')
	ax0.axis([0,50,-mx,mx])
	ax0.legend()

	#fft
	power0,phase0,freq,fr,fi = FFT(t,y[0])
	power1,phase1,freq,fr,fi = FFT(t,y[1])


	#CP
	P0,A0,phi0,Pxyr0,Pxyi0,Freq = CrossPhase(t,y[0],y[1],None,'FFT')
	P1,A1,phi1,Pxyr1,Pxyi1,Freq = CrossPhase(t,y[0],y[1],freq[:-1],'LS')
	
	ax1 = fig.subplot2grid((2,2),(1,0))
	ax1.plot(freq[:-1],power0,color=[1.0,0.0,0.0],label='FFT-0')
	ax1.plot(freq[:-1],power1,color=[1.0,0.5,0.0],label='FFT-1')
	ax1.legend()
	
	ax2 = fig.subplot2grid((2,2),(1,1))
	ax2.plot(freq[:-1],P0,color=[1.0,0.0,0.0],label='CP-FFT')
	ax2.plot(freq[:-1],P1,color=[1.0,0.5,0.0],label='CP-LS')
	ax2.legend()

