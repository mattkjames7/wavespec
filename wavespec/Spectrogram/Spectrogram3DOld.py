import numpy as np
from .Spectrogram import Spectrogram

def Spectrogram3D(t,vx,vy,vz,wind,slip,Freq=None,Method='FFT',
					WindowFunction=None,Param=None,Detrend=True,
					FindGaps=False,GoodData=None,Threshold=0.0,
					Fudge=False,OneSided=True,Tax=None,Steps=None):

	#check that the frequencies exist if we are using LS
	#if Freq is None and Method == 'LS':
	#	print('Please set the Freq keyword before using the LS method')
	#	return



	Nw,LenW,F,xt = Spectrogram(t,vx,wind,slip,Freq,Method,WindowFunction,Param,Detrend,FindGaps,GoodData,Quiet=True,Threshold=Threshold,Fudge=Fudge,OneSided=OneSided,Tax=Tax,Steps=Steps)
	Nw,LenW,F,yt = Spectrogram(t,vy,wind,slip,Freq,Method,WindowFunction,Param,Detrend,FindGaps,GoodData,Quiet=True,Threshold=Threshold,Fudge=Fudge,OneSided=OneSided,Tax=Tax,Steps=Steps)
	Nw,LenW,F,zt = Spectrogram(t,vz,wind,slip,Freq,Method,WindowFunction,Param,Detrend,FindGaps,GoodData,Quiet=True,Threshold=Threshold,Fudge=Fudge,OneSided=OneSided,Tax=Tax,Steps=Steps)
	Nf = F.size - 1
	#need to calculate k vector
	Jxy = xt.Imag*yt.Real - yt.Imag*xt.Real
	Jxz = xt.Imag*zt.Real - zt.Imag*xt.Real
	Jyz = yt.Imag*zt.Real - zt.Imag*yt.Real
	A = np.sqrt(Jxy**2 + Jxz**2 + Jyz**2)	
	kx = Jyz/A
	ky =-Jxz/A
	kz = Jxy/A		
	
	#create an output recarray
	dtype = [('Tspec','float64'),('xSize','int32'),('ySize','int32'),('zSize','int32'),
			('xGood','int32'),('yGood','int32'),('zGood','int32'),
			('xVar','float32'),('yVar','float32'),('zVar','float32'),
			('xPow','float32',(Nf,)),('yPow','float32',(Nf,)),('zPow','float32',(Nf,)),
			('xPha','float32',(Nf,)),('yPha','float32',(Nf,)),('zPha','float32',(Nf,)),
			('xAmp','float32',(Nf,)),('yAmp','float32',(Nf,)),('zAmp','float32',(Nf,)),
			('xReal','float32',(Nf,)),('yReal','float32',(Nf,)),('zReal','float32',(Nf,)),
			('xImag','float32',(Nf,)),('yImag','float32',(Nf,)),('zImag','float32',(Nf,)),
			('kx','float32',(Nf,)),('ky','float32',(Nf,)),('kz','float32',(Nf,))]
	out = np.recarray(Nw,dtype=dtype)
	
	#now fill it up
	names = xt.dtype.names
	for n in names:
		if not n == 'Tspec':
			out['x'+n] = xt[n]
			out['y'+n] = yt[n]
			out['z'+n] = zt[n]

	out.Tspec = xt.Tspec
	out.kx = kx
	out.ky = ky
	out.kz = kz
	
	#clean up previous arrays
	del xt
	del yt
	del zt
	
	return Nw,LenW,F,out
