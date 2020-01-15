import numpy as np

def PolArrow(ax,x0,y0,d=1,r=0.5,color=[0.0,0.0,0.0],linewidth=2.0):
	
	a = d*(45.0 + np.arange(315.0))*np.pi/180.0
	
	x = r*np.sin(a) + x0
	y = r*np.cos(a) + y0
	
	head_width = 0.25 * r
	head_length = 0.25 * r
	
	ax.plot(x,y,color=color,linewidth=linewidth)
	ax.arrow(x[1], y[1], x[0]-x[1], y[0]-y[1], head_width=head_width, head_length=head_length, fc='k', ec='k',color=color)
