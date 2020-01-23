import numpy as np
import DateTimeTools as TT

def UTPlotLabel(fig,axis='x',seconds=False):
	
	if hasattr(fig,'gca'):
		ax=fig.gca()
	else:
		ax = fig
	R = ax.axis()
	mt=ax.xaxis.get_majorticklocs()
	labels=np.zeros(mt.size,dtype='S8')
	for i in range(0,mt.size):
		tmod = mt[i] % 24.0
		hh,mm,ss,ms=TT.DectoHHMM(tmod,True,True,Split=True)

		if seconds:
			utstr=u'{:02n}:{:02n}:{:02n}'.format(hh,mm,ss)
		else:
			if ss >= 30:
				mm+=1
				ss = 0
			if mm > 59:
				hh+=1
				mm=0
			if hh > 23:
				hh = 0
			utstr=u'{:02n}:{:02n}'.format(hh,mm)
		labels[i]=utstr

	labels= np.array(labels).astype('U')
	ax.set_xticks(mt)
	ax.set_xticklabels(labels)
	ax.axis(R)
