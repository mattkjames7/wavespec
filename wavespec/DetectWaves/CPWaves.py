import numpy as np
from .. import Filter
from .. import Fourier

def _BerubeSmooth(Cpha,Cpow,t,f,frange=0.002,trange=1200.0):
	'''
	smooth the spectra using the method defined in Berube et al 2003:
	https://doi.org/10.1029/2002JA009737
	
	
	frange = 0.002 - half the averaging box in Hz
	trange = 1200.0 - half the time averaging box in seconds
	
	'''
	
	#the time input in seconds, so convert trange to seconds
	tr = trange
	fr = frange
	
	if Cpha.shape[0] == t.size:
		nt = t.size
	else:
		nt = t.size - 1
	
	if Cpha.shape[1] == f.size:
		nf = f.size
	else:
		nf = f.size - 1
	
	#output arrays
	pwr_means = np.zeros((nt,nf),dtype='float32') + np.nan
	pha_std = np.zeros((nt,nf),dtype='float32') + np.nan
	pha_means = np.zeros((nt,nf),dtype='float32') + np.nan
	
	for i in range(0,nt):
		uset = np.where((t >= t[i]-tr) & (t <= t[i]+tr))[0]
		t0 = uset.min()
		t1 = uset.max()
		for j in range(0,nf):
			usef = np.where((f >= f[j]-fr) & (f <= f[j]+fr))[0]
			f0 = usef.min()
			f1 = usef.max()			
			
			pwr_means[i,j] = np.nanmean(Cpow[t0:t1+1,f0:f1+1])
			pha_means[i,j] = np.nanmean(Cpha[t0:t1+1,f0:f1+1])
			pha_std[i,j] = np.nanstd(Cpha[t0:t1+1,f0:f1+1],ddof=1)

	return pha_means,pha_std,pha_means


def _SignificantPhase(Cpha):
	
	nt,nf = Cpha.shape
	
	out = np.zeros((nt,nf),dtype='float32') + np.nan
	for i in range(0,nt):
		mean = np.nanmean(Cpha[i])
		std = np.nanstd(Cpha[i],ddof=1)
		for j in range(0,nf):
			if Cpha[i,j] < (mean - std):
				out[i,j] = Cpha[i,j]
				
	return out

def _tTest(Cpha,stds):

	nt,nf = Cpha.shape
	
	ttest = np.zeros((nt,nf),dtype='float32') + np.nan
	Cpha_surv = np.zeros((nt,nf),dtype='float32') + np.nan
	
	for i in range(0,nt):
		for j in range(0,nf):
			tstat = -Cpha[i,j]/stds[i,j]
			ttest[i,j] = tstat
			
			if tstat > 1:
				Cpha_surv[i,j] = Cpha[i,j]
				
	return ttest,Cpha_surv
	
def _Uncertainties(Cpha,F,PR):
	
	nt,nf = Cpha.shape
	
	out = np.zeros((nt,nf),dtype='float32')
	for i in range(0,nt):
		for j in range(0,nf):
			if np.isfinite(Cpha[i,j]):
				
				StartPR = PR[i,j]
				
				if (StartPR > 1) & (j < (nf-2)):
					hiPR = PR[i,j+1]
					cnt = 0
					while (StartPR > 1) & (hiPR > 1) & (j + cnt + 1 < nf):
						StartPR = PR[i,j+cnt]
						hiPR = PR[i,j+cnt+1]
						cnt += 1
						
					if cnt == 0:
						if j > 0:
							out[i,j] = F[j] - F[j-1]
						else:
							out[i,j] = F[j] - F[j+1]
					else:
						x0 = F[j+cnt-1]
						x1 = F[j+cnt]
						y0 = PR[i,j+cnt-1]
						y1 = PR[i,j+cnt]
						m = (y1-y0)/(x1-x0)
						c = y1 - m*x1
						cf = (1-c)/m
						out[i,j] = cf - F[j]
				elif (StartPR < 1) & (i > 0):
					loPR = PR[i,j-1]
					cnt = 0

					while (StartPR < 1) & (loPR < 1) & (j-cnt-1 > 0):
						loPR = PR[i,j-cnt-1]
						StartPR = PR[i,j-cnt]
						cnt += 1
					
					if cnt == 0:
						if j > 0:
							out[i,j] = F[j] - F[j-1]
						else:
							out[i,j] = F[j] - F[j+1]
					else:
						x0 = F[j-cnt]
						x1 = F[j-cnt+1]
						y0 = PR[i,j-cnt]
						y1 = PR[i,j-cnt+1]
						m = (y1-y0)/(x1-x0)
						c = y1 - m*x1
						cf = (1-c)/m
						out[i,j] = F[j] - cf
	return out

def _GetErrBars(t,F,Cpha,unc):
	
	gt,gf = np.where(np.isfinite(Cpha) & (np.abs(unc) < 1))
	
	
	fs = F[gf]
	ts = t[gt]
	us = unc[gt,gf]
	

	return ts,fs,us,gt,gf

def CPWavesFFT(t,Bheq,Bhpol,Window=2400,Slip=300,Highpass=0.00125,**kwargs):
	'''
	This function uses a method based on that used by Wharton et al 2018
	to detect ULF waves using the cross spectra of two magnetometer 
	stations.
	
	
	Inputs
	======
	t : float
		Time axis in seconds
	Bheq : float
		H-component of the magnetic field (nT) measured by the 
		equatorward station
	Bhpol : float
		H-component of the magnetic field (in nT) measured by the 
		poleward station
	Window : float|int
		Window length in seconds
	Slip : float|int
		Gap in seconds between consecutive windows
	Highpass: float
		High-pass filter cutoff (Hz)
		
	Returns
	=======
	
	'''
	
	#filter the data
	res = t[1] - t[0]
	Bhpolf = Filter.Filter(Bhpol,res,high=1/Highpass)
	Bheqf = Filter.Filter(Bheq,res,high=1/Highpass)


	#replace this with CrossSpec
	Nw,F,s_pol = Fourier.Spectrogram(t,Bhpolf,Window,Slip,**kwargs)
	Nw,F,s_eq = Fourier.Spectrogram(t,Bheqf,Window,Slip,**kwargs)
	
	df = F[1:] - F[:-1]
	
	#cross spectrum
	N0 = Window//res
	C = (s_pol.Comp * s_eq.Comp.conj())
	tax = s_pol.Tspec
	tax = np.append(tax,(s_pol.Tspec[-1] + Slip))

	Cpow = np.abs(C) #fudged to be like Sam's
	Cpha = np.angle(C,deg=True)
	
	PR = ((np.abs(s_eq.Comp*N0))**2/(np.abs(s_pol.Comp*N0))**2)

	#get the smoothed spectra
	Cpha_smooth,Cpha_std,Cpow_smooth = _BerubeSmooth(Cpha,Cpow,tax,F)
	
	#largest phase
	Cpha_largest = _SignificantPhase(Cpha_smooth)
	
	#t-test
	ttest,Cpha_surv = _tTest(Cpha_largest,Cpha_std)
	
	#uncertainty
	uncert = _Uncertainties(Cpha_surv,F*1000,PR)
	
	ts,fs,us,ti,fi = _GetErrBars(tax,F,Cpha_surv,uncert)

	#print(ts)
	#print(fs)
	#print(us)
	
	#fill an output dict
	out = {	'Tspec'			:	s_pol.Tspec,
			'Tax'			:	tax,
			'F'				: 	F,
			'C'				:	C,
			'Cpow'			:	Cpow,
			'Cpow_smooth'	:	Cpow_smooth,
			'Cpha'			:	Cpha,
			'Cpha_smooth'	: 	Cpha_smooth,
			'Cpha_largest'	:	Cpha_largest,
			'Cpha_surv'		:	Cpha_surv,
			'Uncert'		:	uncert,
			't'				:	ts,
			'f'				:	fs,
			'ti'			:	ti,
			'fi'			:	fi,
			'u'				:	us/1000.0,
			'ttest'			:	ttest,
			'PR'			:	PR,
			'Tdata'			:	t,
			'Bheq' 			:	Bheqf,
			'Bhpol'			:	Bhpolf}

	return out

def CPWavesFFTSpec(Tspec,F,efft,pfft,N0,**kwargs):
	'''
	This function uses a method based on that used by Wharton et al 2018
	to detect ULF waves using the cross spectra of two magnetometer 
	stations.
	
	
	Inputs
	======
	t : float
		Time axis in seconds
	Bheq : float
		H-component of the magnetic field (nT) measured by the 
		equatorward station
	Bhpol : float
		H-component of the magnetic field (in nT) measured by the 
		poleward station
	Window : float|int
		Window length in seconds
	Slip : float|int
		Gap in seconds between consecutive windows
	Highpass: float
		High-pass filter cutoff (Hz)
		
	Returns
	=======
	
	'''
	

	
	df = F[1:] - F[:-1]
	
	#cross spectrum
	#N0 = Window//res
	C = (pfft* efft.conj())
	tax = Tspec
	Slip = tax[1] - tax[0]
	tax = np.append(tax,(Tspec[-1] + Slip))

	Cpow = np.abs(C) #fudged to be like Sam's
	Cpha = np.angle(C,deg=True)
	
	PR = ((np.abs(efft*N0))**2/(np.abs(pfft*N0))**2)

	#get the smoothed spectra
	Cpha_smooth,Cpha_std,Cpow_smooth = _BerubeSmooth(Cpha,Cpow,tax,F)
	
	#largest phase
	Cpha_largest = _SignificantPhase(Cpha_smooth)
	
	#t-test
	ttest,Cpha_surv = _tTest(Cpha_largest,Cpha_std)
	
	#uncertainty
	uncert = _Uncertainties(Cpha_surv,F*1000,PR)
	
	ts,fs,us,ti,fi = _GetErrBars(tax,F,Cpha_surv,uncert)
	#print(ts)
	#print(fs)
	#print(us)
	
	#fill an output dict
	out = {	'Tspec'			:	Tspec,
			'Tax'			:	tax,
			'F'				: 	F,
			'C'				:	C,
			'Cpow'			:	Cpow,
			'Cpow_smooth'	:	Cpow_smooth,
			'Cpha'			:	Cpha,
			'Cpha_smooth'	: 	Cpha_smooth,
			'Cpha_largest'	:	Cpha_largest,
			'Cpha_surv'		:	Cpha_surv,
			'Uncert'		:	uncert,
			't'				:	ts,
			'f'				:	fs,
			'ti'			:	ti,
			'fi'			:	fi,
			'u'				:	us/1000.0,
			'ttest'			:	ttest,
			'PR'			:	PR,}

	return out

