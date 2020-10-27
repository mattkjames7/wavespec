#include "libfilter.h"


void LanczosKernelLP(int i0, int nt, float *t, float fc, float *k) {
	
	/* start by working out pi*dt/2 */
	int i;
	float t0 = t[i0];
	float *pidt2 = new float[nt];
	for (i=0;i<nt;i++) {
		pidt2[i] = M_PI*fc*(t[i] - t0)/1.5;
	}
	
	/* calculate the function */
	for (i=0;i<nt;i++) {
		k[i] = powf(sinf(pidt2[i])/pidt2[i],2.0)*sinf(3.0*pidt2[i])/(pidt2[i]*3);
	}
	
	/* set the central point to 1.0 because it is NAN */
	k[i0] = 1.0;
	
	/* integrate the function */
	float vol = 0.0;
	for (i=0;i<nt-1;i++) {
		vol += 0.5*(k[i] + k[i+1]);
	}
	
	/* Normalize the function */
	for (i=0;i<nt;i++) {
		k[i]/=vol;
	}
	
	/* clean up */
	delete pidt2;
	
}
	
void LanczosKernelHP(int i0, int nt, float *t, float fc, float *k) {
	
	/* start by working out pi*dt/2 */
	int i;
	float t0 = t[i0];
	float *pidt2 = new float[nt];
	for (i=0;i<nt;i++) {
		pidt2[i] = M_PI*fc*(t[i] - t0)/1.5;
	}
	
	/* calculate the function */
	for (i=0;i<nt;i++) {
		k[i] = powf(sinf(pidt2[i])/pidt2[i],2.0)*sinf(3.0*pidt2[i])/(pidt2[i]*3);
	}
	
	/* set the central point to 1.0 because it is NAN */
	k[i0] = 1.0;
	
	/* integrate the function */
	float vol = 0.0;
	for (i=0;i<nt-1;i++) {
		vol += 0.5*(k[i] + k[i+1]);
	}
	
	/* Normalize the function */
	for (i=0;i<nt;i++) {
		k[i]/=-vol;
	}
	k[i0] = 1.0;
	
	/* clean up */
	delete pidt2;
	
}

float IConv(int nt, float *x, float *k) {
	
	int i;
	float out = 0.0;
	for (i=0;i<nt;i++) {
		out += x[i]*k[i];
	}
	return out;
}

void IFilter(int n, float *t, float *x, float fc, int ftype, float *o) {
	
	/* work out the size of the window (3*cutoff period)
	 * dt is half the window length */
	float dt = 1.5/fc;
	
	/* create some arrays which will say what indices in the time array
	 * fit within the window for each element of the time series */
	int i,j;
	int *i0 = new int[n];
	int *i1 = new int[n];
	float dti;

	/* this algorithm assumes that t is monotonically increasing*/
	for (i=0;i<n;i++) {
		/* start at previous element */
		if (i==0) {
			i0[i] = 0;
			i1[i] = 0;
		} else {
			i0[i] = i0[i-1];
			i1[i] = i1[i-1];
		}
		
		/* loop until within dt of t[i] */
		dti = t[i] - t[i0[i]];
		while ((dti > dt)  && (i0[i] < n-1)) {
			i0[i]++;
			dti = t[i] - t[i0[i]];
		}
		dti = t[i1[i]] - t[i];
		while ((dti <= dt) && (i1[i] < n-1)) {
			i1[i]++;
			dti = t[i1[i]] - t[i];
		}
	}
	
	/* create a temporary array to store the kernel */
	int nt;
	float *k = new float[n];

	/* select the function to use */
	KernelFunc Kernel;
	if (ftype > 0) {
		/* high pass */
		Kernel = &LanczosKernelHP;
	} else {
		/* low pass */
		Kernel = &LanczosKernelLP;
	}
	
	/* loop through each element */
	for (i=0;i<n;i++) {
		/* get the kernel first */
		nt = i1[i] - i0[i] + 1;
		Kernel(i-i0[i],nt,&t[i0[i]],fc,k);
		
		/* now convolve the relevent bits of t and k */
		o[i] = IConv(nt,&x[i0[i]],k);
	}
	
	/* clean up */
	delete i0;
	delete i1;
	delete k;
	
}
