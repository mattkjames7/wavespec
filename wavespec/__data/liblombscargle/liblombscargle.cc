#include "liblombscargle.h"
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
	int i,j;
	double c,s,xc,xs,cc,ss,cs,tau,c_tau,s_tau,c_tau2,s_tau2,cs_tau,w;
	
	for (i=0;i<nf;i++) {
		xc = 0.0;
		xs = 0.0;
		cc = 0.0;
		ss = 0.0;
		cs = 0.0;
		
		w = 2.0*M_PI*f[i];
		
		for (j=0;j<nt;j++) {
			c = cos(w*t[j]);
			s = sin(w*t[j]);
			
			xc += x[j]*c;
			xs += x[j]*s;
			cc += c*c;
			ss += s*s;
			cs += c*s;
		}
		
		//printf("%f %f %f %f %f\n",xc,xs,cc,ss,cs);
		
		tau = atan2(2.0*cs,cc-ss)/(2.0*w);
		c_tau = cos(w*tau);
		s_tau = sin(w*tau);
		c_tau2 = c_tau*c_tau;
		s_tau2 = s_tau*s_tau;
		cs_tau = c_tau*s_tau;
		
		real[i] = xc;
		imag[i] = -xs;
		
		P[i] = 0.5*((pow(c_tau*xc + s_tau*xs,2.0)/(c_tau2*cc + cs_tau*cs + s_tau2*ss)) + (pow(c_tau*xs - s_tau*xc,2.0)/(c_tau2*ss - cs_tau*cs + s_tau2*cc)));
		Pha[i] = atan2(-xs,xc)*180.0/M_PI;
		
	}
	
}

