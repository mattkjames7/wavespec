import numpy as np
from .Fourier.FFT import FFT
from scipy.signal import detrend

def Spectrogram(t,v,wind,slip,WindowFunction=None,Detrend=True,DetectGaps=False,GoodData=None,Quiet=True,LenW=None):



	Tlen = np.size(t)
	if Tlen <= 1:
		return (0,0,0,0,0,0,0,0)

	res = t[1] - t[0]


	if DetectGaps:
		if GoodData is None:
			good = np.zeros(Tlen,dtype='bool')
			gd = np.where(np.isfinite(v))[0]
			good[gd] = True
		else:
			good = GoodData
		st=-1
		ngd=0
		for i in range(0,Tlen):
			if good[i]:
				if st == -1:
					st=i
			else:
				if st != -1:
					if ngd == 0:
						UTi0 = np.array([st])
						UTi1 = np.array([i-1])
					else:
						UTi0 = np.append(UTi0,st)
						UTi1 = np.append(UTi1,i-1)	
					st=-1
					ngd+=1
		if st != -1:
			if ngd == 0:
				UTi0 = np.array([st])
				UTi1 = np.array([i-1])
			else:
				UTi0 = np.append(UTi0,st)
				UTi1 = np.append(UTi1,i-1)	
			st=-1
			ngd+=1
	else:
		ngd = 1
		UTi0 = np.array([0])
		UTi1 = np.array([Tlen-1])

	Nwind=np.zeros(ngd,dtype='int32')
	T0=UTi0.astype('int32')
	Tranges=np.zeros(ngd,dtype='float64')
	for i in range(0,ngd):
		Tranges[i]=t[UTi1[i]]-t[UTi0[i]]
		Nwind[i] = np.int32((Tranges[i]-wind)/slip)+2
	
	
	
	posWind=np.where(Nwind > 0)[0]
	NwTot=np.sum(Nwind[posWind])+ngd-1	
	
	
	
	Tspec=np.zeros(NwTot,dtype='float32')
	if LenW is None:
		LenW=np.int32(np.round(wind/Res))//2
	Pow=np.zeros((NwTot,LenW),dtype='float32')+np.float32(np.nan)
	Pha=np.zeros((NwTot,LenW),dtype='float32')+np.float32(np.nan)
	Real=np.zeros((NwTot,LenW),dtype='float32')+np.float32(np.nan)
	Imag=np.zeros((NwTot,LenW),dtype='float32')+np.float32(np.nan)
	BW=(np.arange(LenW*2+1,dtype='float32')/(LenW*2))/Res
	BW=BW[0:LenW+1]	

	nd=0
	pos=0
	for i in range(0,ngd):
		if nd > 0:
			Tspec[pos]=(Tspec[pos-1]+Tarrays[UTi0[i]]+wind/2.0)/2.0
			pos+=1
		
		if Nwind[i] > 0:
			ng=UTi1[i]-UTi0[i]+1
			good = np.arange(ng)+UTi0[i]
			Tv=v[good]
			TUT=t[good]
			Tax=np.arange(Nwind[i],dtype='float32')*slip + wind/2.0 + TUT[0]
			TaxH=Tax/3600.0
			Tspec[pos:pos+Nwind[i]]=Tax
			for j in range(0,Nwind[i]):
				use0=np.int32(j*slip/Res)
				use=use0+np.arange(LenW*2)
				if np.max(use) >= np.size(Tv):
					bad=1
				else:
					bad=np.where(np.logical_or(Tv[use] < -1e30,np.logical_or( Tv[use] > 1e10 , np.isfinite(Tv[use]) == False)))[0]
				if np.size(bad) == 0:
					if Detrend:
						pwr,pha,b_w,cx,sx=FFT(TUT[use],detrend(Tv[use]),return_complex=True,SmoothComplex=SmoothComplex)
					else:
						pwr,pha,b_w,cx,sx=FFT(TUT[use],Tv[use],return_complex=True,SmoothComplex=SmoothComplex)
					Pow[j+pos]=pwr[0:LenW]
					Pha[j+pos]=pha[0:LenW]
					Real[j+pos]=cx[0:LenW]
					Imag[j+pos]=sx[0:LenW]
				if not Quiet:
					print('\r{:6.2f}%'.format(100.0*np.float32(pos+j+1)/NwTot),end='')
	
			pos+=Nwind[i]
			nd+=1
	if not Quiet:
		print('')
		
	if LowF is None:
		LowF = np.min(BW)
	if HighF is None:
		HighF = np.max(BW)
	
	useF = np.where((BW >= LowF) & (BW <= HighF))[0]
	

	if ReturnComplex:
		return (NwTot,LenW,Tspec,BW[useF],Pow[:,useF[:-1]],Pha[:,useF[:-1]],Real[:,useF[:-1]],Imag[:,useF[:-1]])
	else:
		return (NwTot,LenW,Tspec,BW[useF],Pow[:,useF[:-1]],Pha[:,useF[:-1]])
