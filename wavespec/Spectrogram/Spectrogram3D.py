import numpy as np
from .Spectrogram import Spectrogram

def Spectrogram3D(t,vx,vy,vz,wind,slip,Freq=None,Method='FFT',WindowFunction=None,Param=None,Detrend=True,FindGaps=False,GoodData=None):

	#check that the frequencies exist if we are using LS
	if Freq is None and Method == 'LS':
		print('Please set the Freq keyword before using the LS method')
		return



	Nw,LenW,Freq,xt = Spectrogram(t,vx,wind,slip,Freq,Method,WindowFunction,Param,Detrend,FindGaps,GoodData,Quiet=True)
	Nw,LenW,Freq,yt = Spectrogram(t,vy,wind,slip,Freq,Method,WindowFunction,Param,Detrend,FindGaps,GoodData,Quiet=True)
	Nw,LenW,Freq,zt = Spectrogram(t,vz,wind,slip,Freq,Method,WindowFunction,Param,Detrend,FindGaps,GoodData,Quiet=True)
	
	#need to calculate k vector
	Jxy = xt.Imag*yt.Real - yt.Imag*xt.Real
	Jxz = xt.Imag*zt.Real - zt.Imag*xt.Real
	Jyz = yt.Imag*zt.Real - zt.Imag*yt.Real
	A = np.sqrt(Jxy**2 + Jxz**2 + Jyz**2)	
	kx = Jyz/A
	ky =-Jxz/A
	kz = Jxy/A		
	
	#create an output recarray
	dtype = [('Tspec','float32'),
			('xPow','float32',(LenW,)),('yPow','float32',(LenW,)),('zPow','float32',(LenW,)),
			('xPha','float32',(LenW,)),('yPha','float32',(LenW,)),('zPha','float32',(LenW,)),
			('xAmp','float32',(LenW,)),('yAmp','float32',(LenW,)),('zAmp','float32',(LenW,)),
			('xReal','float32',(LenW,)),('yReal','float32',(LenW,)),('zReal','float32',(LenW,)),
			('xImag','float32',(LenW,)),('yImag','float32',(LenW,)),('zImag','float32',(LenW,)),
			('kx','float32',(LenW,)),('ky','float32',(LenW,)),('kz','float32',(LenW,))]
	out = np.recarray(Nw,dtype=dtype)
	
	#now fill it up
	names = xt.dtype.names
	for n in names:
		out['x'+n] = xt[n]
		out['y'+n] = yt[n]
		out['z'+n] = zt[n]
	out.kx = kx
	out.ky = ky
	out.kz = kz
	
	#clean up previous arrays
	del xt
	del yt
	del zt
	
	return Nw,LenW,Freq,out
