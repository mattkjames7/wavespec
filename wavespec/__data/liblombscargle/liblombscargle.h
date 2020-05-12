#ifndef __liblombscargle_h__
#define __liblombscargle_h__
#include <math.h>
#include <stdio.h>
#include "ComplexLS.h"

#endif
using namespace std;

extern "C" {
	void LombScargle(double *t, double *x, int n, double *f, int nf, double *P, double *A, double *phi, double *a, double *b, double Threshold, bool Fudge);
	void LombScargleSam(float *t, float *x, int n, float *f, int nf, float *P, float *A, float *phi, float *a, float *b);
}
