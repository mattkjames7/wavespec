import numpy as np
import matplotlib.pyplot as plt
from ..Fourier import FFT
from ..LombScargle import LombScargle

def TestLS(A=[1.0,2.0],f=[0.04,0.1],phi=[0.0,90.0],Backend='C++'):
	'''
	Test the Lomb-Scargle routine against the FFT.
	 
	'''
	#create the time series
	n = np.size(A)
	t = np.arange(0,51.0)
	yi = np.zeros((n,51),dtype='float32')
	for i in range(0,n):
		yi[i] = A[i]*np.cos(2*np.pi*f[i]*t + np.pi*phi[i]/180.0)
	y = np.sum(yi,axis=0)
	
	mx = np.max(np.abs(y))
	
	fig = plt
	fig.figure(figsize=(8,8))
	ax0 = fig.subplot2grid((2,1),(0,0))
	for i in range(0,n):
		ax0.plot(t,yi[i],label='$A_{'+'{:d}'.format(i)+'}$='+'{:3.1f}'.format(A[i])+' $f_{'+'{:d}'.format(i)+'}$='+'{:4.2f}'.format(f[i])+r' $\phi_{'+'{:d}'.format(i)+'}$='+'{:5.1f}'.format(phi[i])) 
	ax0.plot(t,y,label='Combined')
	ax0.set_xlabel('$t$')
	ax0.set_ylabel('$y$')
	ax0.plot([0.0,50.0],[0.0,0.0],color=[0.0,0.0,0.0],linestyle=':')
	ax0.axis([0,50,-mx,mx])
	ax0.legend()

	#fft
	power,phase,freq,fr,fi = FFT(t,y)

	fpk = np.where((power[1:-1] > power[:-2]) & (power[1:-1] > power[2:]))[0] + 1

	#LS
	P,A,phi,a,b = LombScargle(t,y,freq[1:-1],Backend)
	
	lpk = np.where((P[1:-1] > P[:-2]) & (P[1:-1] > P[2:]))[0] + 1
	lfreq = freq[1:]

	ax1 = fig.subplot2grid((2,2),(1,0))
	ax1.plot(freq[:-1],power,color=[1.0,0.0,0.0])
	
	nfpk = fpk.size
	for i in range(0,nfpk//2):
		ax1.text(0.1,0.9-i*0.1,'$A_{'+'{:d}'.format(i)+'}$='+'{:3.1f}'.format(np.sqrt(power[fpk[i]])*2)+' $f_{'+'{:d}'.format(i)+'}$='+'{:4.2f}'.format(freq[fpk[i]])+r' $\phi_{'+'{:d}'.format(i)+'}$='+'{:5.1f}'.format(phase[fpk[i]]*180.0/np.pi))
	
	ax2 = fig.subplot2grid((2,2),(1,1))
#	ax2.plot(freq[1:-1],(A/2)**2,color=[1.0,0.0,0.0])
	ax2.plot(freq[1:-1],P,color=[1.0,0.0,0.0])
	
	nlpk = lpk.size
	for i in range(0,nlpk//2):
		ax2.text(0.1,0.9-i*0.1,'$A_{'+'{:d}'.format(i)+'}$='+'{:3.1f}'.format(np.sqrt(P[lpk[i]])*2)+' $f_{'+'{:d}'.format(i)+'}$='+'{:4.2f}'.format(lfreq[lpk[i]])+r' $\phi{'+'{:d}'.format(i)+'}$='+'{:5.1f}'.format(phi[lpk[i]]*180.0/np.pi))
	
	
