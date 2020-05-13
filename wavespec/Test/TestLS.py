import numpy as np
import matplotlib.pyplot as plt
from ..Fourier import FFT
from ..LombScargle import LombScargle

def TestLS(A=[1.0,2.0],f=[0.04,0.1],phi=[0.0,90.0],Backend='C++',Threshold=0.01):
	'''
	Test the Lomb-Scargle routine against the FFT.
	 
	'''
	#create the time series
	n = np.size(A)
	t = np.arange(0,50.0)
	yi = np.zeros((n,50),dtype='float32')
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
	ax0.plot([0.0,49.0],[0.0,0.0],color=[0.0,0.0,0.0],linestyle=':')
	ax0.axis([0,49,-mx,mx])
	ax0.legend()

	#fft
	power,Af,phase,fr,fi,freq = FFT(t,y,OneSided=True,Threshold=Threshold)

	fpk = np.where((power[1:-1] > power[:-2]) & (power[1:-1] > power[2:]))[0] + 1

	#LS
	P,A,phi,a,b = LombScargle(t,y,freq[:-1],Backend,Threshold=Threshold)
	
	lpk = np.where((P[1:-1] > P[:-2]) & (P[1:-1] > P[2:]))[0] + 1
	print(freq)

	ax1 = fig.subplot2grid((2,2),(1,0))
	#ax1.plot(freq[:-1],power,color=[1.0,0.0,0.0])
	ax1.stem(freq[:-1],power)
	
	nfpk = fpk.size
	for i in range(0,nfpk):
		ax1.text(0.02,0.93-i*0.06,'$A_{'+'{:d}'.format(i)+'}$='+'{:3.1f}'.format(np.sqrt(power[fpk[i]])*2)+' $f_{'+'{:d}'.format(i)+'}$='+'{:4.2f}'.format(freq[fpk[i]])+r' $\phi_{'+'{:d}'.format(i)+'}$='+'{:5.1f}'.format(phase[fpk[i]]*180.0/np.pi),transform=ax1.transAxes)
	
	ax2 = fig.subplot2grid((2,2),(1,1))
#	ax2.plot(freq[1:-1],(A/2)**2,color=[1.0,0.0,0.0])
	#ax2.plot(freq[1:-1],P,color=[1.0,0.0,0.0])
	ax2.stem(freq[:-1],P)
	
	nlpk = lpk.size
	for i in range(0,nlpk):
		ax2.text(0.02,0.93-i*0.06,'$A_{'+'{:d}'.format(i)+'}$='+'{:3.1f}'.format(np.sqrt(P[lpk[i]])*2)+' $f_{'+'{:d}'.format(i)+'}$='+'{:4.2f}'.format(freq[lpk[i]])+r' $\phi{'+'{:d}'.format(i)+'}$='+'{:5.1f}'.format(phi[lpk[i]]*180.0/np.pi),transform=ax2.transAxes)	
	
	
	ax0.set_title('Time Series')
	ax0.set_xlabel('Time (s)')
	ax0.set_ylabel('$y$')
	
	
	ax1.set_title('FFT')
	ax1.set_xlabel('Frequency (Hz)')
	ax1.set_ylabel('Power')
	
	
	ax2.set_title('L-S')
	ax2.set_xlabel('Frequency (Hz)')
	ax2.set_ylabel('Power')
	
	
	return fr,fi,a,b



def TestLS2(Amp=0.0,Phase=0.0,DC=0.0,BE=['Python','C++','Sam'],OneSided=True):
	'''
	Do a quick test to see whether FFT and L-S match up. Note that if
	doing a two-sided FFT, then LS amplitude should be twice the FFT
	amplitude and the LS power should be 4 times the FFT power.
	
	'''
	
	#create a time series
	t = np.arange(32.0)
	x = DC + Amp*np.cos(2.0*np.pi*0.25*t + Phase) + 0.1*np.random.rand(t.size)
	
	
	
	
	#FFT
	P0,A0,p0,r0,i0,f = FFT(t,x,OneSided=OneSided)
	f = f[:-1]
	

	fig = plt
	fig.figure(figsize=(8,11))
	
	#Power
	ax0 = fig.subplot2grid((5,1),(0,0))
	ax0.plot(f,P0,label='FFT')

	ax0.set_title('Power')
	
	#Amp
	ax1 = fig.subplot2grid((5,1),(1,0))
	ax1.plot(f,A0,label='FFT')

	ax1.set_title('Amplitude')
	
	#Phase
	ax2 = fig.subplot2grid((5,1),(2,0))
	ax2.plot(f,p0,label='FFT')

	ax2.set_title('Phase')
	
	#Real
	ax3 = fig.subplot2grid((5,1),(3,0))
	ax3.plot(f,r0,label='FFT')

	ax3.set_title('Real')
	
	#Imag
	ax4 = fig.subplot2grid((5,1),(4,0))
	ax4.plot(f,i0,label='FFT')

	ax4.set_title('Imaginary')

	
	for i in range(0,len(BE)):
		P1,A1,p1,r1,i1 = LombScargle(t,x,f,BE[i],Fudge=True)
		ax0.plot(f,P1,label=BE[i])
		ax1.plot(f,A1,label=BE[i])
		ax2.plot(f,p1,label=BE[i])
		ax3.plot(f,r1,label=BE[i])
		ax4.plot(f,i1,label=BE[i])
		
	ax0.legend()
	ax1.legend()
	ax2.legend()
	ax3.legend()
	ax4.legend()

	fig.tight_layout()
