#ifndef _fint1_h
#define _fint1_h

#include <rsf.h>

typedef struct Vint1 *vint1;
typedef struct Fint1 *fint1;

float* fint1_coeff (fint1 fnt, int n);
float* vint1_coeff (vint1 fnt, int n, int dim);
vint1 vint1_init (int nw, int n1, int dim);
fint1 fint1_init (int nw, int n1);
void vint1_set (vint1 fnt, float** dat);
void fint1_set (fint1 fnt, float* dat);
void fint1_close (fint1 fnt);
void vint1_close (vint1 fnt);
float fint1_apply (fint1 fnt, int i, float x, bool der);
void vint1_apply (vint1 fnt, int i, float x, bool der, float* f);
void stretch(fint1 fnt, 
	     float (*map)(float),
	     int n1, float d1, float o1,
	     int n2, float d2, float o2,
	     float *trace);

#endif

/* 	$Id: fint1.h,v 1.3 2004/04/01 02:12:42 fomels Exp $	 */
