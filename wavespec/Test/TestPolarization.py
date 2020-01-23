import numpy as np
import matplotlib.pyplot as plt
from ..Tools.ArrowLine import ArrowLine
from ..Tools.Polarization2D import Polarization2D

def TestPolarization(xPow=2.0,xPhase=0.0,yPow=1.0,yPhase=40.0):
	
	xPhase = xPhase*np.pi/180.0
	yPhase = yPhase*np.pi/180.0
	
	
	Ax,Ay,Axi,Aeta,psi,e,direction = Polarization2D(xPow,xPhase,yPow,yPhase)
	
	
	t = np.arange(361.0)*np.pi*2
	x = Ax*np.cos(t/360.0 + xPhase)
	y = Ay*np.cos(t/360.0 + yPhase)
	
	xi = Axi*np.cos(t/360.0)
	eta = Aeta*np.cos(t/360.0 + np.pi/2.0)
	
	fig = plt
	fig.figure(figsize=(8,11))
	ax0 = fig.subplot2grid((3,1),(0,0))
	ax1 = fig.subplot2grid((3,1),(1,0))
	ax2 = fig.subplot2grid((3,2),(2,0))
	ax3 = fig.subplot2grid((3,2),(2,1))

	ax0.plot(t,x,color=[1.0,0.5,0.0],label='$x$')
	ax0.plot(t,y,color=[1.0,0.0,0.5],label='$y$')
	ax0.legend()
	
	ax1.plot(t,xi,color=[1.0,0.5,0.5],label=r'$\xi$')
	ax1.plot(t,eta,color=[0.0,0.5,0.5],label=r'$\eta$')
	ax1.legend()
	
	
	x_xi = Axi*np.cos(psi)
	y_xi = Axi*np.sin(psi)
	x_eta = Aeta*np.cos(psi+np.pi/2.0)
	y_eta = Aeta*np.sin(psi+np.pi/2.0)
	
	
	
	ax2.set_aspect(1.0)
	ArrowLine(ax2,x,y,color=[0.0,0.0,0.0],HeadWidth=0.1,HeadLength=0.2)
	ax2.plot([-Ax,Ax,Ax,-Ax,-Ax],[-Ay,-Ay,Ay,Ay,-Ay],color=[0.0,0.0,0.0])
	ax2.plot([-x_xi,x_xi],[-y_xi,y_xi],color=[1.0,0.5,0.5],label=r'$\xi$')
	ax2.plot([-x_eta,x_eta],[-y_eta,y_eta],color=[0.0,0.5,0.5],label=r'$\eta$')
	ax2.set_aspect(1.0)
	R = ax2.axis()
	ax2.axis(R)
	ax2.plot([R[0],R[1]],[0.0,0.0],color=[0.0,0.0,0.0],linestyle=':')
	ax2.plot([0.0,0.0],[R[2],R[3]],color=[0.0,0.0,0.0],linestyle=':')
	ax2.legend(loc='upper left')
	ax2.set_xlabel('$x$')
	ax2.set_ylabel('$y$')

	if direction == 1:
		dirs = 'RH'
	elif direction == -1:
		dirs = 'LH'
	else:
		dirs = 'Linear'
	
	ax3.axis([0,1,0,1])
	ax3.axis('off')
	ax3.text(0.1,0.9,r'$\psi$ = {:6.2f}'.format(psi*180.0/np.pi))
	ax3.text(0.1,0.8,r'$\delta$ = {:6.2f}'.format((yPhase-xPhase)*180.0/np.pi))
	ax3.text(0.1,0.7,r'$e$ = {:6.4f}'.format(e))
	ax3.text(0.1,0.6,r'Handedness: {:s}'.format(dirs))
	ax3.text(0.1,0.5,r'$A_x$ = {:4.2f}'.format(Ax))
	ax3.text(0.1,0.4,r'$A_y$ = {:4.2f}'.format(Ay))	
	ax3.text(0.1,0.3,r'$A_{\xi}$ = ' + '{:4.2f}'.format(Axi))
	ax3.text(0.1,0.2,r'$A_{\eta}$ = ' + '{:4.2f}'.format(Aeta))
