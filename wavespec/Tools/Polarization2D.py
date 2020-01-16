import numpy as np

def Polarization2D(xPow,xPha,yPow,yPha):
	'''
	Calculates approximate wave amplitudes, eccentricity and polarization
	direction using the power and phase in two dimensions.
	
	
	The calculation of the polarization ellipse comes from section 1.4.2
	of Born and Wolf: "Principles of Optics" 
	
	Inputs
	======
	xPow : Wave power in x direction 
	xPha : Wave phase in x direction
	yPow : Wave power in y direction 
	yPha : Wave phase in y direction
	
	Returns
	=======
	Ax : Amplitude in x direction
	Ay : Amplitude in y direction
	Axi : Amplitude in xi direction
	Aeta : Amplitude in eta direction
	e : eccentricity
	direction : +1 = right-handed, -1 = left-handed, 0 = perfectly linear
	
	NOTE: The handedness here is defined by the direction in which the 
	wave rotates about the z-axis, where right-handedness (left-handedness)
	corresponds to the wave rotating anti-clockwise (clockwise) about 
	the z-axis when viewed by an observer at +z looking down towards the
	x-y plane.
	
	'''
	
	
	#first we need to know the difference in phase
	delta = yPha - xPha
	
	#also the amplitudes
	Ax = 2*np.sqrt(xPow)
	Ay = 2*np.sqrt(yPow)
	
	#the auxiliary angle
	alpha = np.arctan2(Ay,Ax)
	
	#the angle between the x/y axes and the axes of the ellipse
	psi = np.arctan(np.tan(2.0*alpha)*np.cos(delta))/2.0
	
	#calculate the axes of the ellipse
	cospsi = np.cos(psi)
	sinpsi = np.sin(psi)
	cos2psi = cospsi**2
	sin2psi = sinpsi**2
	Const = 2*Ax*Ay*cospsi*sinpsi*np.cos(delta)
	
	Axi = np.sqrt(cos2psi*Ax**2 + sin2psi*Ay**2 + Const)
	Aeta = np.sqrt(sin2psi*Ax**2 + cos2psi*Ay**2 - Const)
	
	#calculate eccentricity
	e = np.sqrt(1 - (Aeta**2)/(Axi**2))
	
	
	#hopefully this bit will correctly define the direction (handedness) 
	#of the wave around the z axis
	direction = np.sign(((delta % (2*np.pi)))-np.pi)
	
	
	return (Ax,Ay,Axi,Aeta,psi,e,direction)	
