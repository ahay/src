#include <math.h>

#include <rsf.h>

#include "pick0.h"
#include "eno.h"

static int n1, n2;
static eno *ent;
static float pmax, pmin;

void pick0_init (int n1_in, int n2_in, int order)
{
    int i2;

    n1 = n1_in;
    n2 = n2_in;

    ent = (eno*) sf_alloc (n2,sizeof(eno));
    for (i2 = 0; i2 < n2; i2++) {
	ent[i2] = eno_init (order, n1);
    }

    pmin = 2.*n1;
    pmax = - 2.*n1;
}

void pick0_set (int i2, float* dip)
{
    int i1;

    eno_set (ent[i2], dip);
    for (i1 = 0; i1 < n1; i1++) {
	if      (pmin > dip[i1]) pmin = dip[i1];
	else if (pmax < dip[i1]) pmax = dip[i1];
    }
}

void pick0_close (void)
{
    int i2;
    
    for (i2 = 0; i2 < n2; i2++) {
	eno_close (ent[i2]);
    }
    free (ent);
}

void pick0_step (float t0, float* t)
{
    int i, i2;
    float k1, k2;

    t[0] = t0;
    for (i2 = 0; i2 < n2-1; i2++) {
	i = floor(t0); eno_apply (ent[i2], i, t0-i, &k1, NULL, FUNC);
	if (k1 < pmin) k1 = pmin;
	if (k1 > pmax) k1 = pmax;
	t0 += k1; 
	i = floor(t0); eno_apply (ent[i2+1], i, t0-i, &k2, NULL, FUNC);
	if (k2 < pmin) k2 = pmin;
	if (k2 > pmax) k2 = pmax;
	t0 += 0.5*(k2-k1);
	t[i2+1] = t0;
    }
}

void pick0_step0 (float t0, float* t)
{
    int i, i2;
    float k1, k2;

    t[0] = t0;
    for (i2 = 0; i2 < n2-1; i2++) {
	i = floor(t0); eno_apply (ent[i2], i, t0-i, &k1, NULL, FUNC);
	if (k1 < pmin) k1 = pmin;
	if (k1 > pmax) k1 = pmax;
	t0 += k1;
	t[i2+1] = t0;
    }
}

void pick0_delta (int m2, float* t)
{
    int i, i2, i1;
    float ti, k1, k2;
    
    for (i1=0; i1 < n1; i1++) {
	ti=i1;
	for (i2=m2; i2 > 0; i2--) {
	    i = floor(ti); eno_apply (ent[i2], i, ti-i, &k1, NULL, FUNC);
	    if (k1 < pmin) k1 = pmin;
	    if (k1 > pmax) k1 = pmax;
	    ti -= k1; 
	    i = floor(ti); eno_apply (ent[i2-1], i, ti-i, &k2, NULL, FUNC);
	    if (k2 < pmin) k2 = pmin;
	    if (k2 > pmax) k2 = pmax;
	    ti += 0.5*(k1-k2);
	}
	t[i1]=ti-i1;
    }
}
