#ifndef __liblombscargle_h__
#define __liblombscargle_h__
#include <math.h>
#include <stdio.h>
#include "ComplexLS.h"

#endif
using namespace std;

extern "C" {
	
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
						double *b, double Threshold, bool Fudge);
						

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
	void LombScargleSam(double *t, double *x, int n, double *f, int nf, 
					double *P, double *A, double *phi, double *a, double *b);
}
