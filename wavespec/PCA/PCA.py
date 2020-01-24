import numpy as np
from scipy.signal import detrend


def GetPrincipalComponents(x,y,z):

	mux = np.mean(x)
	muy = np.mean(y)
	muz = np.mean(z)
	
	sxx = np.sum((x-mux)**2)
	sxy = np.sum((x-mux)*(y-muy))
	sxz = np.sum((x-mux)*(z-muz))
	syy = np.sum((y-muy)**2)
	syz = np.sum((y-muy)*(z-muz))
	szz = np.sum((z-muz)**2)
	
	
	data = np.array([x-mux,y-muy,z-muz]).transpose()
	
	#ScatterMatrix=np.array([[sxx,sxy,sxz],[sxy,syy,syz],[sxz,syz,szz]])
	#CoVarMatrix=np.cov([x-mux,y-muy,z-muz])
	#Russell calls the scatter matrix a variance matrix
	#the covariance matrix = variance/(N-1)
	#according to a website these should both provide the same eigenvalues
	#though in this case they don't seem to
	ScatterMatrix = np.array([[sxx,sxy,sxz],[sxy,syy,syz],[sxz,syz,szz]])
	CoVarMatrix = np.cov(data.T)
	
	# eigenvectors and eigenvalues for the from the covariance matrix
	EigValCov, EigVecCov = np.linalg.eig(CoVarMatrix)
	
	srt = (EigValCov.argsort())[::-1]
	eigs = EigValCov[srt]
	vecs = EigVecCov[:,srt]
	
	return (eigs,vecs.T)
	
def RotateToPrincipalComponents(x,y,z):
	
	eigs,vecs = GetPrincipalComponents(x,y,z)
	
	nx = x*vecs[0,0] + y*vecs[0,1] + z*vecs[0,2]
	ny = x*vecs[1,0] + y*vecs[1,1] + z*vecs[1,2]
	nz = x*vecs[2,0] + y*vecs[2,1] + z*vecs[2,2]

	return (nx,ny,nz)


def SlidingPCA(t,x,y,z,Wind,Slip,Detrend=True):
	
	Tlen = np.size(t)

	Trange = t[Tlen-1]-t[0]
	Nwind = np.int32((Trange-Wind)/Slip)+1

	MinVec = np.zeros((Nwind,3),dtype='float32')+np.float32(np.nan)
	MidVec = np.zeros((Nwind,3),dtype='float32') +np.float32(np.nan)
	MaxVec = np.zeros((Nwind,3),dtype='float32') +np.float32(np.nan)
	Eigs = np.zeros((Nwind,3),dtype='float32')+np.float32(np.nan)
	
	Twind = np.arange(Nwind,dtype='float32')*Slip + Wind/2.0 + t[0]
	
	LenW = np.int32(Wind/0.05)/2
	for i in range(0,Nwind):
		use0 = np.int32(i*Slip/0.05)
		use = use0+np.arange(LenW*2)
		if np.max(use) >= np.size(t):
			bad = 1
		else:
			bad = np.where((np.isfinite(x) == False) | (np.isfinite(y) == False) | (np.isfinite(z) == False))[0]
		if np.size(bad) == 0:
			if Detrend:
				eigs,vecs = GetPrincipalComponents(detrend(x[use]),detrend(y[use]),detrend(z[use]))
			else:
				eigs,vecs = GetPrincipalComponents(x[use],y[use],z[use])	
			Eigs[i] = eigs
			MaxVec[i] = vecs[0]
			MidVec[i] = vecs[1]
			MinVec[i] = vecs[2]
			
	return (Twind,Eigs,MinVec,MidVec,MaxVec)



def CovarianceMatrix(x):
	'''
	Input x must have shape (m,n), where m is the number of samples and n is the number of features
	
	Output is a (n,n) matrix
	'''
	mu = np.mean(x,axis=0)
	cov = np.dot((x-mu).T,(x-mu))/(x.shape[0]-1)
	
	return cov
	

def EigenValues(x):
	
	cov = CovarianceMatrix(x)
	EigVal,EigVec = np.linalg.eig(cov)
	
	srt = np.argsort(EigVal)[::-1]
	
	return EigVal[srt],EigVec[:,srt]
	

def PCA(x,nComp=None,EigVal=None,EigVec=None,ReturnEigVal=False,ReturnEigVec=False):
	if EigVal is None or EigVec is None:
		EigVal,EigVec = EigenValues(x)
	ev = np.abs(EigVal)
	
	if not nComp is None:
		#either a real or integer 
		if np.int32(nComp) == nComp:
			#easy - just select number of components
			EigVal1 = EigVal[:np.int32(nComp)]
			EigVec1 = EigVec[:,:np.int32(nComp)]
			ev1 = ev[:np.int32(nComp)]
		else:
			#hopefully a floating point, in which case assume it is a number between 0.0 and 1.0
			cs = np.cumsum(ev)/np.sum(ev)
			use = np.where(cs >= nComp)[0]
			if use.size > 1:
				ind = use[1]
			elif use.size == 1:
				ind = use[0]
			else:
				ind = ev.size-1
			EigVal1 = EigVal[:ind+1]
			EigVec1 = EigVec[:,:ind+1]
			ev1 = ev[:ind+1]
	else:
		EigVec1 = EigVec
	x_transformed = np.dot(x,EigVec1)
	
	if not ReturnEigVec and not ReturnEigVal:
		out = x_transformed
	else:
		out = [x_transformed]
	if ReturnEigVal:
		out.append(EigVal)
		
	if ReturnEigVec:
		out.append(EigVec)
		
	return out
			
