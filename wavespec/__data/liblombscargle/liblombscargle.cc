#include "liblombscargle.h"

/* Calculates angular frequency */
/***********************************************************************
 * NAME : _angfreq(f,nf,w)
 * 
 * DESCRIPTION : 
 * 		Calculate the angular frequency from the frequency.
 * 
 * INPUTS :
 * 		double *f		Array of frequencies (Hz).
 * 		int nf			Length of frequency array.
 * 
 * RETURNS : 
 * 		double *w		Angular frequency array (rad s^-1).
 * 
 **********************************************************************/
void _angfreq(double *f, int nf, double *w) {
	int i;
	for (i=0;i<nf;i++) {
		w[i] = 2*M_PI*f[i];
	}
} 


/* Calculate variance in data array, x, assuming mean has already been
 * subtracted */
 /***********************************************************************
 * NAME : _variance(x,n)
 * 
 * DESCRIPTION : 
 * 		Calculate the variance in an array.
 * 
 * INPUTS :
 * 		double *x			value array - mean
 * 		int n				number of elements 		
 * 	
 * RETURNS : 
 * 		double var			Variance.
 * 
 **********************************************************************/
double _variance(double *x, int n) {
	int i;
	double o2 = 0.0;
	for (i=0;i<n;i++) {
		o2 += x[i]*x[i];
	}
	return o2/(n-1);
}

/* Calculates Tau for a given frequency */
/***********************************************************************
 * NAME : _tau(t,n,w)
 * 
 * DESCRIPTION : 
 * 		Calculate the tau parameter for a frequency.
 * 
 * INPUTS :
 * 		double *t 		Time array
 * 		int n			Length of array
 * 		double w		Angular frequency
 * 
 * RETURNS : 
 * 		double tau
 * 
 * 
 **********************************************************************/
double _tau(double *t, int n, double w) {
	int i;
	double wt;
	double ss2w = 0.0;
	double sc2w = 0.0;
	for (i=0;i<n;i++) {
		wt = 2*w*t[i];
		ss2w += sin(wt);
		sc2w += cos(wt);
	}
	return atan2(ss2w,sc2w)/(2*w);
}

/* Calculates w(t_i - tau) */
/***********************************************************************
 * NAME : _wtT(t,n,w,tau,wtT)
 * 
 * DESCRIPTION : 
 * 		Calculates w(t_i - tau) 
 * 
 * INPUTS :
 * 		double *t		Time array.
 * 		int n			number of elements
 * 		double w		angular frequency
 * 		double tau		tau
 * 
 * OUTPUTS : 
 * 		double *wtT		
 * 
 **********************************************************************/
void _wtT(double *t, int n, double w, double tau, double *wtT) {
	int i;
	for (i=0;i<n;i++) {
		wtT[i] = w*(t[i] - tau);
		//wtT[i] = w*(t[i]);
	}
}

/* Calculates cos(w(t_i - tau)) */
/***********************************************************************
 * NAME : _coswtT(n,wtT,coswtT)
 * 
 * DESCRIPTION : 
 * 		Calculates cos(w(t_i - tau))
 * 
 * INPUTS :
 * 		int	n		Number of elements.
 * 		double *wtT			wtT
 * 
 * OUTPUTS : 
 * 		double *coswtT
 * 
 **********************************************************************/
void _coswtT(int n, double *wtT, double *coswtT) {
	int i;
	for (i=0;i<n;i++) {
		coswtT[i] = cos(wtT[i]);
	}
}
/* Calculates sin(w(t_i - tau)) */
/***********************************************************************
 * NAME : _coswtT(n,wtT,coswtT)
 * 
 * DESCRIPTION : 
 * 		Calculates sin(w(t_i - tau))
 * 
 * INPUTS :
 * 		int	n		Number of elements.
 * 		double *wtT			wtT
 * 
 * OUTPUTS : 
 * 		double *sinwtT
 * 
 **********************************************************************/
void _sinwtT(int n, double *wtT, double *sinwtT) {
	int i;
	for (i=0;i<n;i++) {
		sinwtT[i] = sin(wtT[i]);
	}
}


 
/* Calculates the summations used for working out a and b*/
/***********************************************************************
 * NAME : _sums(x,n,coswtT,sinwtT,syc,sys,sc2,ss2)
 * 
 * DESCRIPTION : 
 * 		Calcualte the sums used.
 * 
 * INPUTS :
 * 		double *x		Data array
 * 		int n			Number of elements
 * 		double *coswtT	cos(w(t_i - tau))
 * 		double *sinwtT	sin(w(t_i - tau))
 * 
 * OUTPUTS : 
 * 		double * syc	Sum of x*coswtT
 * 		double * sys	Sum of x*sinwtT
 * 		double * sc2	Sum of coswtT*coswtT
 * 		double * ss2	Sum of sinwtT*sinwtT
 * 		
 * 
 **********************************************************************/
void _sums(double *x, int n, double *coswtT, double *sinwtT, double *syc, double *sys, double *sc2, double *ss2) {
	int i;
	syc[0] = 0.0;
	sys[0] = 0.0;
	sc2[0] = 0.0;
	ss2[0] = 0.0;
	for (i=0;i<n;i++) {
		syc[0] += x[i]*coswtT[i];
		sys[0] += x[i]*sinwtT[i];
		sc2[0] += coswtT[i]*coswtT[i];
		ss2[0] += sinwtT[i]*sinwtT[i];
	}
}

/***********************************************************************
 * NAME : _mean(x,n)
 * 
 * DESCRIPTION : 
 *  	Calculate the mean of an array.
 * 
 * INPUTS :
 * 		double *x		Array to be averaged.
 * 		int n			Length of the array.
 * 
 * RETURNS : 
 * 		double mean		Mean of array x.
 **********************************************************************/
double _mean(double *x, int n) {
	int i;
	double tmp = 0.0;
	for (i=0;i<n;i++) {
		tmp += x[i];
	}
	return tmp/((double) n);
}


/***********************************************************************
 * NAME : LombScargle(t,x,n,f,nf,P,A,phi,a,b,Threshold,Fudge)
 * 
 * DESCRIPTION : 
 *  	Basic Lomb Scargle analysis with additional phase calculation.
 * 
 * Inputs:
 * 		double *t			time array, length n
 * 		double *x			value array, length array
 * 		int n				number of time elements
 *		double *f			frequency array (Hz), length nf
 * 		int nf 				number of frequencies
 * 
 * Outputs:	
 * 		double *P			Periodogram output, length nf
 * 		double *A			Amplitude output, length nf
 * 		double *Pha			Wave phase in radians, length nf
 * 		double *a			Parameter a, length nf
 * 		double *b			Parameter b, length nf
 * 		double Threshold	Threshold amplitude level, below which all
 * 							LS data will be set to zero in order to 
 * 							reduce noise. If Threshold=0.0 then nothing
 * 							happens
 * 		bool Fudge			This will enable/disable a fudge for when 
 * 							there is likely to be an almost 0/0 
 * 							division, but not quite due to floating 
 * 							point errors, relatively untested but it 
 * 							provides the correct	values at half the 
 * 							nyquist frequency, matching FFT.
 */
void LombScargle(	double *t, double *x, int n, double *f, int nf, 
					double *P, double *A, double *phi, double *a, 
					double *b, double Threshold, bool Fudge) {
	/*create some variables*/
	int i,j;
	double *w = new double[nf];
	double *wtT = new double[n];
	double *coswtT = new double[n];
	double *sinwtT = new double[n];
	double o2;
	double tau;
	double syc, sys, sc2, ss2;
	double rt2n = sqrt(2.0/((float) n));
	double a2b2;
	
	/*Convert frequency (Hz) to angular frequency, omega (rad/s)*/
	_angfreq(f,nf,w);
	
	/* Calculate variance in data set */
	//o2 = _variance(x,n);
	
	/* Loop through each frequency */
	for (i=0;i<nf;i++) {
		if (f[i] == 0.0) {
			A[i] = _mean(x,n);
			P[i] = A[i]*A[i];
			phi[i] = 0.0;
			a[i] = A[i];
			b[i] = 0.0;
		} else {
			/* Get tau and then w*(t[i] - tau)*/
			tau = _tau(t,n,w[i]);
			_wtT(t,n,w[i],tau,wtT);
		
			/* Calculate sine and cosine of w*(t[i] - tau) */
			_sinwtT(n,wtT,sinwtT);
			_coswtT(n,wtT,coswtT);
		
			/* calculate sums */
			_sums(x,n,coswtT,sinwtT,&syc,&sys,&sc2,&ss2);
			
			/* a and b */
			a[i] = 0.5*rt2n*syc/sqrt(sc2);
			b[i] = 0.5*rt2n*sys/sqrt(ss2);
			if (Fudge) {
				/* set to zero if there is likely to be a problem with
				 * the fractions above because of small numbers*/
				if ((fabs(syc) < 1e-8) && (fabs(sc2) < 1e-8)) {
					a[i] = 0.0;
				}
				if ((fabs(sys) < 1e-8) && (fabs(ss2) < 1e-8)) {
					b[i] = 0.0;
				}
			}
			
			/* Power */
			P[i] = (a[i]*a[i] + b[i]*b[i])*4.0;
			
			/* Amplitude */
			A[i] = sqrt(P[i]);
			
			/* Phase */
			phi[i] = fmod(((-atan2(b[i],a[i]) - w[i]*tau) + 3*M_PI),(2.0*M_PI)) - M_PI;
			
			/* Real and imaginary terms*/
			a[i] = A[i]*cos(phi[i])/2.0;
			b[i] = A[i]*sin(phi[i])/2.0;
			
			
		}
		
		
	}
	
	/* apply a threshold to stuff */
	if (Threshold > 0.0) {
		for (i=0;i<nf;i++) {
			if (A[i] < Threshold) {
				A[i] = 0.0;
				P[i] = 0.0;
				phi[i] = 0.0;
				a[i] = 0.0;
				b[i] = 0.0;
			}				
		}
	}
	
	/* Deallocate arrays */
	delete w;
	delete wtT;
	delete sinwtT;
	delete coswtT;
	
}

/***********************************************************************
 * NAME : LombScargleSam(t,x,n,f,nf,P,A,phi,a,b)
 * 
 * DESCRIPTION : 	This function is a wrapper for Samuel Wharton's
 * 		complex Lomb Scargle code which essentially calculates a DFT for 
 * 		irregularly spaced data.
 * 
 * INPUTS : 
 * 		double *t			Time array (s).
 * 		double *x			Data to be transformed.
 * 		int n				Length of the time array.
 * 		double *f			Desired frequencies (Hz).
 * 		int nf				Length of the omega array.
 * 
 * OUTPUTS : 
 * 		double *P			Output power.
 * 		double *A			Output amplitudes.
 * 		double *phi			Output phase.
 * 		double *a 			Parameter a
 * 		double *b			Parameter b
 * 
 * Parameters a and b are related to the real and imaginary components
 * of the FFT.
 * 
 * ********************************************************************/
void LombScargleSam(	double *t, double *x, int n, double *f, int nf, 
			double *P, double *A, double *phi, double *a, double *b) {

	/*convert frequency to omega*/
	int i,j;
	double *w = new double[nf];
	_angfreq(f,nf,w);

	/*call Sam's code*/
	ComplexLS(n,nf,t,x,w,A,phi,a,b);

	/*calculate power and the real/imaginary components*/
	for (i=0;i<nf;i++) {
		P[i] = A[i]*A[i];
		a[i] = A[i]*cos(phi[i])/2.0;
		b[i] = A[i]*sin(phi[i])/2.0;
	}
	
	/*deallocate omega*/
	delete w;

}
