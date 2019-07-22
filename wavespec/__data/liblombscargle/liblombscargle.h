#ifndef __liblombscargle_h__
#define __liblombscargle_h__
#include <math.h>
#include <stdio.h>
#endif
using namespace std;

extern "C" {
	void LombScargle(double *t, double *x, int nt, double *f, int nf, double *P, double *Pha, double *real, double *imag);
}
