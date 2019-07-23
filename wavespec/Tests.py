import numpy as np
import matplotlib.pyplot as plt
from .LombScargle.LombScargle import LombScargle
import scipy.signal as signal

def TestLS():
	'''
	Perform a quick test on the Lomb Scargle routine, compare to scipy. 
	
	'''
	
	#create the time series (code from scipy docs)
	A0 = 2.0
	f0 = np.random.rand()*0.4+0.1
	w0 = 2*np.pi*f0
	phi = 2*np.pi*np.random.rand()
	nin = 1000
	nout = 100000
	frac_points = 0.9
	r = np.random.rand(nin)
	x = np.linspace(0.01,10*np.pi,nin)
	x = x[r >= frac_points]
	y = A0*np.sin(w0*x + phi)
	w = np.linspace(0.01,10,nout)
	f = w/(2.0*np.pi)
	
	P,A,Phi,a,b = LombScargle(x,y,f,'C++')
	Pp,Ap,Phip,ap,bp = LombScargle(x,y,f,'Python')
	Pnorm = P/P.max()
	Ppnorm = Pp/Pp.max()
	
	pgram = signal.lombscargle(x,y,w,normalize=True)


	fig = plt 
	fig.figure()
	#plot power (normalized)
	ax0 = fig.subplot2grid((2,1),(0,0))
	ax0.plot(f,Pnorm,color=[0.0,1.0,0.0],label='C++')
	ax0.plot(f,Ppnorm,color=[0.0,0.0,1.0],linestyle='--',label='Python')
	ax0.plot(f,pgram,color=[1.0,0.0,0.0],linestyle=':',label='scipy')
	ax0.legend()
	
	#find peak and plot wave
	xwave = np.linspace(x.min(),x.max(),1000)
	pk = np.where(P == np.nanmax(P))[0][0]
	fwave = f[pk]
	phiwave = Phi[pk]
	awave = A[pk]
	ywave = awave*np.cos(2*np.pi*fwave*xwave + phiwave)
	print(fwave,phiwave)
	pkp = np.where(Pp == np.nanmax(Pp))[0][0]
	fpwave = f[pkp]
	phipwave = Phip[pk]
	apwave = Ap[pk]
	
	ypwave = apwave*np.cos(2*np.pi*fpwave*xwave + phipwave)
	print(fpwave,phipwave)
	
	ax1 = fig.subplot2grid((2,1),(1,0))
	ax1.scatter(x,y,color=[0.0,0.0,0.0,0.2],label='Data')
	ax1.plot(xwave,ywave,color=[0.0,1.0,0.0],label='C++')
	ax1.plot(xwave,ypwave,color=[0.0,0.0,1.0],linestyle='--',label='Python')
	ax1.legend()
