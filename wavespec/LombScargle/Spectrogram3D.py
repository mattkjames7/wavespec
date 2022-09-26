import numpy as np
from .Spectrogram import Spectrogram
from ..Tools.kVector import kVector

def Spectrogram3D(t,vx,vy,vz,wind,slip,**kwargs):
	'''
	Calculate 3D spectrogram of sets of time series data on a common
	time axis. Alternatively, spatial 3D spectrogram of 3 datasets with
	a common spatial axis (just swap s and Hz for m and m^-1).
	
	Inputs
	======
	t : float
		Time axis (s).
	vx : float
		x-component data.
	vy : float
		y-component data.
	vz : float
		z-component data.
	wind : float
		Window length (s).
	slip : float
		Time from one window start to the next (s).
	
	Keyword Args
	============
	CombineComps : bool
		If True, some additional Pow, Pha, Amp and Comp fields included 
		in the output by combining the spectra from each pair of 
		components.
	For other keywords see the Spectrogram() function.
		
	Returns
	=======
	Nw : int
		Number of windows.
	Freq : float
		Frequency array (Hz) - one element longer than the number of
		elements in each spectrum.
	out : numpy.recarray
		All of the outputs from the Spectrogram() function with the 
		component name as a prefix (e.g. 'Pow' -> 'xPow'). Also, all
		three components of the wave unit vector (k) are included.
		
	
	
	'''
	CombineComps = kwargs.get('CombineComps',False)

	#Calculate the three sets of spectra
	Nw,F,xt = Spectrogram(t,vx,wind,slip,**kwargs)
	Nw,F,yt = Spectrogram(t,vy,wind,slip,**kwargs)
	Nw,F,zt = Spectrogram(t,vz,wind,slip,**kwargs)
	Nf = F.size -1
	
	#need to calculate k vector
	kx,ky,kz = kVector(xt.Comp,yt.Comp,zt.Comp)
	
	#create an output recarray data type
	dtype = [	('Tspec','float64'),
				('xSize','int32'),
				('ySize','int32'),
				('zSize','int32'),
				('xGood','int32'),
				('yGood','int32'),
				('zGood','int32'),
				('xVar','float32'),
				('yVar','float32'),
				('zVar','float32'),
				('xPow','float32',(Nf,)),
				('yPow','float32',(Nf,)),
				('zPow','float32',(Nf,)),
				('xPha','float32',(Nf,)),
				('yPha','float32',(Nf,)),
				('zPha','float32',(Nf,)),
				('xAmp','float32',(Nf,)),
				('yAmp','float32',(Nf,)),
				('zAmp','float32',(Nf,)),
				('xComp','complex64',(Nf,)),
				('yComp','complex64',(Nf,)),
				('zComp','complex64',(Nf,)),
				('kx','float32',(Nf,)),
				('ky','float32',(Nf,)),
				('kz','float32',(Nf,))]
			
	#combine some components
	if CombineComps:
		dtypec = [	('xyComp','complex64',(Nf,)),
					('yzComp','complex64',(Nf,)),
					('zxComp','complex64',(Nf,)),
					('xyPow','float32',(Nf,)),
					('yzPow','float32',(Nf,)),
					('zxPow','float32',(Nf,)),
					('xyPha','float32',(Nf,)),
					('yzPha','float32',(Nf,)),
					('zxPha','float32',(Nf,)),
					('xyAmp','float32',(Nf,)),
					('yzAmp','float32',(Nf,)),
					('zxAmp','float32',(Nf,)),]
		for dc in dtypec:
			dtype.append(dc)
			
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
	
	if CombineComps:
		cc = ['xy','yz','zx']
		for c in cc:
			c0 = c[0]
			c1 = c[1]
			Comp = out[c0+'Comp'] * np.conjugate(out[c1+'Comp'])
			out[c+'Comp'] = Comp
			out[c+'Amp'] = np.abs(Comp)
			out[c+'Pow'] = out[c+'Amp']**2
			out[c+'Pha'] = np.arctan2(Comp.imag,Comp.real)
	
	#clean up previous arrays
	del xt
	del yt
	del zt
	
	return Nw,F,out
