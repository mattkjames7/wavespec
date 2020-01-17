#include "liblombscargle.h"

/* Calculates angular frequency */
void _angfreq(double *f, int nf, double *w) {
	int i;
	for (i=0;i<nf;i++) {
		w[i] = 2*M_PI*f[i];
	}
} 

/* Calculate variance in data array, x, assuming mean has already been
 * subtracted */
double _variance(double *x, int n) {
	int i;
	double o2 = 0.0;
	for (i=0;i<n;i++) {
		o2 += x[i]*x[i];
	}
	return o2/(n-1);
}

/* Calculates Tau for a given frequency */
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
void _wtT(double *t, int n, double w, double tau, double *wtT) {
	int i;
	for (i=0;i<n;i++) {
		//wtT[i] = w*(t[i] - tau);
		wtT[i] = w*(t[i]);
	}
}

/* Calculates cos(w(t_i - tau)) */
void _coswtT(int n, double *wtT, double *coswtT) {
	int i;
	for (i=0;i<n;i++) {
		coswtT[i] = cos(wtT[i]);
	}
}
/* Calculates sin(w(t_i - tau)) */
void _sinwtT(int n, double *wtT, double *sinwtT) {
	int i;
	for (i=0;i<n;i++) {
		sinwtT[i] = sin(wtT[i]);
	}
}


 
/* Calculates the summations used for working out a and b*/
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

double _mean(double *x, int n) {
	int i;
	double tmp = 0.0;
	for (i=0;i<n;i++) {
		tmp += x[i];
	}
	return tmp/((double) n);
}


/* Basic Lomb Scargle analysis with additional phase calculation
 * 
 * Inputs:
 * t		time array, length nt
 * x		value array, length array
 * n		number of time elements
 * f		frequency array (Hz), length nf
 * nf 		number of frequencies
 * 
 * Outputs:
 * P		Periodogram output, length nf
 * A		Amplitude output, length nf
 * Pha		Wave phase in radians, length nf
 * a		Parameter a, length nf
 * b		Parameter b, length nf
 */
void LombScargle(double *t, double *x, int n, double *f, int nf, double *P, double *A, double *phi, double *a, double *b) {
	/*create some variables*/
	int i,j;
	double *w = new double[nf];
	double *wtT = new double[n];
	double *coswtT = new double[n];
	double *sinwtT = new double[n];
	double o2;
	double tau;
	double syc, sys, sc2, ss2;
	double rt2n = sqrt(2.0/n);
	double a2b2;
	
	/*Convert frequency (Hz) to angular frequency, omega (rad/s)*/
	_angfreq(f,nf,w);
	
	/* Calculate variance in data set */
	o2 = _variance(x,n);
	
	/* Loop through each frequency */
	for (i=0;i<nf;i++) {
		if (f[i] == 0.0) {
			A[i] = _mean(x,n);
			phi[i] = NAN;
			a[i] = NAN;
			b[i] = NAN;
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
			a[i] = rt2n*syc/sqrt(sc2);
			b[i] = rt2n*sys/sqrt(ss2);
			
			/* Periodogram */
			a2b2 = a[i]*a[i] + b[i]*b[i];
			//P[i] = a2b2*((double) n)/(4.0*o2);
			P[i] = a2b2/4;
			
			/* Amplitude */
			A[i] = sqrt(a2b2);
			
			/* Phase */
			phi[i] = -atan2(b[i],a[i]);
		}
	}
	
	/* Deallocate arrays */
	delete w;
	delete wtT;
	delete sinwtT;
	delete coswtT;
	
}

