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


def PlotPolEllipse(xPow,xPhase,yPow,yPhase,fig=None,maps=[1,1,0,0]):


	xPhase = xPhase*np.pi/180.0
	yPhase = yPhase*np.pi/180.0
	
	
	Ax,Ay,Axi,Aeta,psi,e,direction = Polarization2D(xPow,xPhase,yPow,yPhase)
	
	if direction == 1:
		dirs = 'RH'
	elif direction == -1:
		dirs = 'LH'
	else:
		dirs = 'Linear'
	
	
	t = np.arange(361.0)*np.pi*2
	x = Ax*np.cos(t/360.0 + xPhase)
	y = Ay*np.cos(t/360.0 + yPhase)
	
	xi = Axi*np.cos(t/360.0)
	eta = Aeta*np.cos(t/360.0 + np.pi/2.0)

	x_xi = Axi*np.cos(psi)
	y_xi = Axi*np.sin(psi)
	x_eta = Aeta*np.cos(psi+np.pi/2.0)
	y_eta = Aeta*np.sin(psi+np.pi/2.0)
	
		
	
	if fig is None:
		fig = plt
		fig.figure()
	if hasattr(fig,'Axes'):	
		ax = fig.subplot2grid((maps[1],maps[0]),(maps[3],maps[2]))
	else:
		ax = fig

	ax.set_aspect(1.0)
	ArrowLine(ax,x,y,color=[0.0,0.0,0.0],HeadWidth=0.1,HeadLength=0.2)
	ax.plot([-Ax,Ax,Ax,-Ax,-Ax],[-Ay,-Ay,Ay,Ay,-Ay],color=[0.0,0.0,0.0])
	ax.plot([-x_xi,x_xi],[-y_xi,y_xi],color=[1.0,0.5,0.5],label=r'$\xi$')
	ax.plot([-x_eta,x_eta],[-y_eta,y_eta],color=[0.0,0.5,0.5],label=r'$\eta$')
	ax.set_aspect(1.0)
	R = ax.axis()
	Rm = np.max(np.abs(np.array(R)))
	fs = 5.0
	ax.axis([-Rm,Rm,-Rm,Rm])
	ax.plot([R[0],R[1]],[0.0,0.0],color=[0.0,0.0,0.0],linestyle=':')
	ax.plot([0.0,0.0],[R[2],R[3]],color=[0.0,0.0,0.0],linestyle=':')
	ax.legend(loc='upper right',fontsize=fs)
	ax.set_xlabel('$x$')
	ax.set_ylabel('$y$')
	
	ax.text(0.05,0.875,r'$\psi$ = {:6.2f}'.format(psi*180.0/np.pi),transform=ax.transAxes,fontsize=fs)
	ax.text(0.05,0.75,r'$\delta$ = {:6.2f}'.format((yPhase-xPhase)*180.0/np.pi),transform=ax.transAxes,fontsize=fs)
	ax.text(0.05,0.625,r'$e$ = {:6.4f}'.format(e),transform=ax.transAxes,fontsize=fs)
	ax.text(0.05,0.5,r'Handedness: {:s}'.format(dirs),transform=ax.transAxes,fontsize=fs)
	ax.text(0.05,0.375,r'$A_x$ = {:4.2f}'.format(Ax),transform=ax.transAxes,fontsize=fs)
	ax.text(0.05,0.25,r'$A_y$ = {:4.2f}'.format(Ay),transform=ax.transAxes,fontsize=fs)	
	ax.text(0.05,0.125,r'$A_{\xi}$ = ' + '{:4.2f}'.format(Axi),transform=ax.transAxes,fontsize=fs)
	ax.text(0.05,0.0,r'$A_{\eta}$ = ' + '{:4.2f}'.format(Aeta),transform=ax.transAxes,fontsize=fs)

	
	
	return ax

def _RemovexTicks(ax):
	
	xt = ax.get_xticks()
	xtl = len(xt)*['']
	ax.set_xticks(xt)
	ax.set_xticklabels(xtl)

def _RemoveyTicks(ax):
	
	xt = ax.get_yticks()
	xtl = len(xt)*['']
	ax.set_yticks(xt)
	ax.set_yticklabels(xtl)

def TestMultiPolarization():
	
	#list of amplitudes and phases
	a = 0.707
	xA = np.array([1.0,a,0.0,a,1.0,1.0])
	yA = np.array([0.0,a,1.0,a,1.0,1.0])
	xP = np.array([0.0,0.0,0.0,0.0,90.0,45.0])
	yP = np.array([0.0,0.0,0.0,180.0,0.0,0.0])
	nw = 6
	
	#make a long time series
	t = np.arange(nw*1000).astype('float64')
	f = 0.01
	x = np.zeros(nw*1000,dtype='float64')
	y = np.zeros(nw*1000,dtype='float64')
	z = np.zeros(nw*1000,dtype='float64')
	for i in range(0,nw):
		ind = np.arange(1000) + 1000*i
		x[ind] = xA[i]*np.cos(2*np.pi*f*t[ind] + xP[i]*np.pi/180.0)
		y[ind] = yA[i]*np.cos(2*np.pi*f*t[ind] + yP[i]*np.pi/180.0)
	
	#spectral analysis
	from ..Spectrogram.Spectrogram3D import Spectrogram3D
	Nw,F,Spec = Spectrogram3D(t,x,y,z,400,50,CombineComps=True)
	
		
	#create figure
	fig = plt.figure(figsize=(8,11))
	
	
	#plot polarization ellipses

	ax0 = []
	for i in range(0,nw):
		axt = PlotPolEllipse(xA[i]**2,xP[i],yA[i]**2,yP[i],fig=plt,maps=[nw,8,i,0])
		_RemovexTicks(axt)
		_RemoveyTicks(axt)
		ax0.append(axt)
		
			
	#plot time series
	ax1 = plt.subplot2grid((8,1),(1,0))
	ax1.plot(t,x,color='red',label='$x$')
	ax1.plot(t,y,color='green',label='$y$')
	ax1.legend()
	ax1.set_xlim(t[0],t[-1])
	ax1.set_xlabel('$t$')
	_RemovexTicks(ax1)

	
	#plot x power
	from ..Spectrogram.SpectrogramPlotter import SpectrogramPlotter
	ax2 = SpectrogramPlotter(Spec.Tspec,F*1000,Spec.xPow,fig=plt,maps=[1,8,0,2],zlabel='$x$ Power')
	ax2.set_ylim(0,f*2000)
	ax2.set_xlim(t[0],t[-1])
	_RemovexTicks(ax2)
	
	#plot y power
	ax3 = SpectrogramPlotter(Spec.Tspec,F*1000,Spec.yPow,fig=plt,maps=[1,8,0,3],zlabel='$y$ Power')
	ax3.set_ylim(0,f*2000)
	ax3.set_xlim(t[0],t[-1])	
	_RemovexTicks(ax3)
	
	#plot combined power (with peaks)
	ax4 = SpectrogramPlotter(Spec.Tspec,F*1000,Spec.xyPow,fig=plt,maps=[1,8,0,4],zlabel='$x \times y$ Power')
	ax4.set_ylim(0,f*2000)
	ax4.set_xlim(t[0],t[-1])
	_RemovexTicks(ax4)
	
	#plot combined power (with peaks)
	ax5 = SpectrogramPlotter(Spec.Tspec,F*1000,Spec.xPow+Spec.yPow,fig=plt,maps=[1,8,0,5],zlabel='$x + y$ Power')
	ax5.set_ylim(0,f*2000)
	ax5.set_xlim(t[0],t[-1])
	_RemovexTicks(ax5)
	
	fi = np.array([p.argmax() for p in (Spec.xPow+Spec.yPow)])
	ts = Spec.Tspec
	ax5.scatter(ts,F[fi]*1000.0,marker='.',color='lime')
	
	xPow = Spec.xPow[:,fi]
	xPha = Spec.xPha[:,fi]
	yPow = Spec.yPow[:,fi]
	yPha = Spec.yPha[:,fi]
	
	Ax,Ay,Axi,Aeta,psi,e,direction = Polarization2D(xPow,xPha,yPow,yPha)
	
	#plot frequency, psi, e
	ax6 = plt.subplot2grid((8,1),(6,0))
	ax6.plot(ts,psi,color='red',label=r'$\psi')
	ax6.plot(ts,e,color='green',label='$e$')
	ax6.set_xlim(t[0],t[-1])
	ax6.set_ylim([-np.pi,np.pi])
	ax6.hlines([-np.pi/2,0.0,np.pi/2.0],t[0],t[-1],color='black',linestyle=':')
	ax6r = ax6.twinx()
	ax6r.plot(ts,F[fi]*1000.0,color='orange',label='$f$')
	ax6r.set_ylim(0.0,f*2000)
#	ax6.legend()
	_RemovexTicks(ax6)
	
	#plot Ax,Ay and what they should be
	ax7 = plt.subplot2grid((8,1),(7,0))
	ax7.plot(ts,Ax,color='red',label='$A_x$')
	ax7.plot(ts,Ay,color='green',label='$A_y$')
	_RemovexTicks(ax7)
	#plot Axi,Aeta and what they should be


	plt.subplots_adjust(top=0.97,bottom=0.1,left=0.05,right=0.95,hspace=0.0)
