import numpy as np
from .. import Fourier
from .. import LombScargle

def Spectrogram3D(t,vx,vy,vz,wind,slip,**kwargs):

	Method = kwargs.get('Method','FFT')

	if Method == 'FFT':
		return Fourier.Spectrogram3D(t,vx,vy,vz,wind,slip,**kwargs)
	elif Method == 'LS':
		return LombScargle.Spectrogram3D(t,vx,vy,vz,wind,slip,**kwargs)
	else:
		print('Method not supported')
		return None
