import numpy as np
import matplotlib.pyplot as plt

def ArrowLine(fig,X,Y,color=[0.0,0.0,0.0],zorder=None,Spacing=20,SingleArrow=False,EndArrow=False,linewidth=1,linestyle='-',HeadWidth=0.05, HeadLength=0.1,Reverse=False,**kwargs):
	
	if Reverse == True:
		x=X[::-1]
		y=Y[::-1]
	else:
		x=X
		y=Y
	

	
	if type(fig) == type(plt):
		ax=fig.gca()
	else:
		#assume that axes have been provided instead of a figure
		ax = fig
		
	fig.plot(x,y,color=color,zorder=zorder,linewidth=linewidth,linestyle=linestyle,**kwargs)

	if EndArrow:
		n=1
		i0=np.array([np.size(x)])-2
		i1=i0+1
	elif SingleArrow:
		n=1
		i0=np.array([np.size(x)//2])-1
		i1=i0+1
	else:
		l=np.size(x)
		n=l//Spacing -1
		i0=(np.arange(n)+1)*Spacing+np.int32(0.5*Spacing)
		i1=i0+1
			
	for i in range(0,n):
		ax.arrow(x[i0[i]], y[i0[i]], x[i1[i]]-x[i0[i]], y[i1[i]]-y[i0[i]], head_width=HeadWidth, head_length=HeadLength,color=color,zorder=zorder)
