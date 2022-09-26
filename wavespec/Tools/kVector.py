import numpy as np


def kVector(Fx,Fy,Fz):
	'''
	Calculate the wave vector k using 3D Fourier spectra.
	
	Inputs
	======
	Fx : complex
		x-component FFT/LS spectra/spectrum.
	Fy : complex
		y-component FFT/LS spectra/spectrum.
	Fz : complex
		z-component FFT/LS spectra/spectrum.
		
	Returns
	=======
	kx : float
		x-component of the wave unit vector
	ky : float
		y-component of the wave unit vector
	kz : float
		z-component of the wave unit vector
		
	
	'''
	
	Jxy = Fx.imag*Fy.real - Fy.imag*Fx.real
	Jxz = Fx.imag*Fz.real - Fz.imag*Fx.real
	Jyz = Fy.imag*Fz.real - Fz.imag*Fy.real
	A = np.sqrt(Jxy*Jxy + Jxz*Jxz + Jyz*Jyz)	
	kx = Jyz/A
	ky =-Jxz/A
	kz = Jxy/A		
	
	return kx,ky,kz
	
