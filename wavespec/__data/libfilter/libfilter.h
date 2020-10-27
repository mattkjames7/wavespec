#ifndef __LIBFILTER_H__
#define __LIBFILTER_H__
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#endif
using namespace std;

typedef void (*KernelFunc)(int,int,float*,float,float*);

extern "C" {

	void LanczosKernelLP(int i0, int nt, float *t, float fc, float *k);
	void LanczosKernelHP(int i0, int nt, float *t, float fc, float *k);
	float IConv(int nt, float *x, float *k);
	void IFilter(int n, float *t, float *x, float fc, int ftype, float *o);
}
