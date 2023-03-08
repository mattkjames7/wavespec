import numpy as np
from .. import Globals
from ..Filter.Filter import Filter
from .. import Fourier
from .. import CrossSpec

def _ReadData():
    fname = Globals.ModulePath + 'cpdata.bin'
    return RT.ReadRecarray(fname)

def _FilteredData(Highpass=0.00125):

    data = _ReadData()
    res = 10.0
    data.muo = Filter(data.muo,res,high=1/Highpass)
    data.pel = Filter(data.pel,res,high=1/Highpass)

    return data


def _SamsXSpec(t,Be,Bp,Window,Slip):

	#replace this with CrossSpec
	Nw,F,s_pol = Fourier.Spectrogram(t,Bp,Window,Slip,**kwargs)
	Nw,F,s_eq = Fourier.Spectrogram(t,Be,Window,Slip,**kwargs)
	
	df = F[1:] - F[:-1]
	
	#cross spectrum
	N0 = Window//res
	C = (s_pol.Comp * s_eq.Comp.conj())
	tax = s_pol.Tspec
	tax = np.append(tax,(s_pol.Tspec[-1] + Slip))

	Cpow = np.abs(C) #fudged to be like Sam's
	Cpha = np.angle(C,deg=True)
        
	return tax,F,C,Cpow,Cpha

def _MyXSpec(t,Be,Bp,Window,Slip):
    


def TestCPEigenfreqs():
