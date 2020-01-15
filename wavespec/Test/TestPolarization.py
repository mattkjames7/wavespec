import numpy as np
import matplotlib.pyplot as plt
from ..Tools.ArrowLine import ArrowLine
from ..Tools.Polarization2D import Polarization2D

def TestPolarization(xPow=1.0,xPhase=0.0,yPow=1.0,yPhase=90.0):
	
	xPhase = xPhase*np.pi/180.0
	yPhase = yPhase*np.pi/180.0
	
	
	Ax,Ay,Axi,Aeta,psi,e,direction = Polarization2D(xPow,xPhase,xPow,xPhase)
	
	t = np.arange(361.0)*np.pi*2
	x = Ax*np.cos(t/360.0 + xPhase)
	y = Ay*np.cos(t/360.0 + yPhase)
	
	xi = Axi*np.cos(t/360.0)
	eta = Aeta*np.cos(t/360.0 + np.pi/2.0)
	
	fig = plt
	fig.figure(figsize=(8,11))
	ax0 = fig.subplot2grid((3,1),(0,0))
	ax1 = fig.subplot2grid((3,1),(1,0))
	ax2 = fig.subplot2grid((3,1),(2,0))

	ax0.plot(t,x,color=[1.0,0.5,0.0],label='$x$')
	ax0.plot(t,y,color=[1.0,0.0,0.5],label='$y$')
	ax0.legend()
	
	ax1.plot(t,xi,color=[1.0,0.5,0.5],label=r'$\xi$')
	ax1.plot(t,eta,color=[0.0,0.5,0.5],label=r'$\eta$')
	ax1.legend()
	
	ax2.set_aspect(1.0)
	ArrowLine(ax2,x,y,color=[0.0,0.0,0.0])
	ax2.plot([-Ax,Ax,Ax,-Ax],[-Ay,-Ay,Ay,Ay],color=[0.0,0.0,0.0])
	ax2.plot([0.0,0.0],[xi*np.cos(psi),xi*np.sin(psi)],color=[1.0,0.5,0.5])
	ax2.plot([0.0,0.0],[-xi*np.cos(psi),-xi*np.sin(psi)],color=[1.0,0.5,0.5])
	ax2.plot([0.0,0.0],[eta*np.cos(psi+np.pi/2.0),eta*np.sin(psi+np.pi/2.0)],color=[0.0,0.5,0.5])
	ax2.plot([0.0,0.0],[-eta*np.cos(psi+np.pi/2.0),-eta*np.sin(psi+np.pi/2.0)],color=[0.0,0.5,0.5])
	ax2.set_xlabel('$x$')
	ax2.set_ylabel('$y$')
	
