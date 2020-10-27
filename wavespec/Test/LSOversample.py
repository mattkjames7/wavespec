import numpy as np
import matplotlib.pyplot as plt
from ..LombScargle.LombScargle import LombScargle
from ..Tools.WindowFunctions import ApplyWindowFunction,WindowScaleFactor

def LSOversample(f=0.005,NoiseLevel=1.0,
				Backend='C++',WindowFunction=None,Param=None,
				flim=None):
	'''
	Test the ability of the code to create an oversampled l-s 
	periodogram.
	
	Inputs
	======
	f : float
		Frequency in Hz
	NoiseLevel : float
		Amount of noise to add to data
	'''
	
	#create the time series, 1s time interval
	n = 1000
	t = np.arange(n)
	x = np.cos(2*np.pi*f*t) + NoiseLevel*(2*np.random.rand(n)-1.0)
	o = np.ones(n)
	
	#work out the frequencies
	freq0 = (np.arange(n+1,dtype='float32')/(np.float32(n)))[:n//2]
	freq1 = np.linspace(freq0[0],freq0[-1],n*2)
	
	
	#now do both L-S periodograms
	P0,_,_,_,_ = LombScargle(t,x,freq0,Backend,WindowFunction,Param)
	P1,_,_,_,_ = LombScargle(t,x,freq1,Backend,WindowFunction,Param)
	
	
	#create the plot
	fig = plt
	fig.figure(figsize=(6,6))
	ax0 = fig.subplot2grid((2,1),(0,0))
	ax1 = fig.subplot2grid((2,1),(1,0))
	
	if WindowFunction is None:
		wflab = 'No window'
	else:
		wflab = 'WF: '+WindowFunction
		
	sf = WindowScaleFactor(WindowFunction,Param)
	
	ax0.plot(t,x,color='black')
	ax0.plot(t,ApplyWindowFunction(t,x,WindowFunction,Param),color='red')
	ax0.plot(t,ApplyWindowFunction(t,o,WindowFunction,Param),color='orange',label=wflab)
	ax0.plot(t,ApplyWindowFunction(t,-o,WindowFunction,Param),color='orange')
	
	ax0.set_xlabel('Time (s)')
	ax0.set_title('$f$ = {:6.2f} mHz'.format(1000.0*f))
	ax0.legend()

	ax1.plot(freq0,P0/(sf**2),label='FFT Frequencies',color='blue')
	ax1.plot(freq1,P1/(sf**2),label='4x Oversampled',color='orange')
	if flim is None:
		flim = [0.0,f*2]
	ax1.set_xlim(flim)
	ax1.legend()
	
	ax1.set_xlabel('Frequency (Hz)')
	ax1.set_ylabel('Power')
	
	ylim = ax1.get_ylim()
	ax1.plot([f,f],ylim,color='red',linestyle='--')
	ax1.set_ylim(ylim)
	ax1.text(0.02,0.95,'Backend: '+Backend,transform=ax1.transAxes,ha='left',va='center')

	fig.tight_layout()
