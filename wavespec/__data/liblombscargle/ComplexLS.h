#ifndef __COMPLEXLS_H__
#define __COMPLEXLS_H__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#endif
/***********************************************************************
 * NAME : ComplexLS(lent,lenw,time,data,omega,amplitude,phase,a,b)
 * 
 * DESCRIPTION : 	This function calculates the complex Lomb-Scargle
 * 		periodogram - essentially a DFT for irregularly spaced data.
 * 		This specific function was originally written by Samuel Wharton
 * 		then adapted by me...
 * 
 * INPUTS : 
 * 		int lent			Length of the time array.
 * 		int lenw			Length of the omega array.
 * 		double *time			Time array.
 * 		double *data			Data to be transformed.
 * 		double *omega		Desired angular frequencies.
 * 
 * OUTPUTS : 
 * 		double *amplitude	Output amplitudes.
 * 		double *phase		Output phase.
 * 		double *a 			Parameter a
 * 		double *b			Parameter b
 * 
 * Parameters a and b are related to the real and imaginary components
 * of the FFT.
 * 
 * ********************************************************************/
void ComplexLS(	int lent, int lenw, double *time, double *data, 
				double *omega, double *amplitude, double *phase, 
				double *a, double *b);
