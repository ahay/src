#include <float.h>
#include <math.h>

#include <rsf.h>

#include "monof.h"

/* destroys data */
float monof(float *data, int i0, int niter, float a0, int nk, 
	    float dk, bool verb)
{
    int ik, iter;
    float f2, fe, ee, lfe, le2, l2e2, k, da, aa, f, e2;
    float le, eps, num, den, r2, a, e;
    
    eps = 10.*FLT_EPSILON;
    eps *= eps;

    f2 = 0.;
    for (ik=0; ik < nk; ik++) {
	f = data[ik];
	f2 += f*f;
    }

    a = a0;
    aa = 1.;
	
    for (iter = 0; iter < niter; iter++) {
	ee = eps;
	fe = lfe = le2 = l2e2 = 0.;
	for (ik = 0; ik < nk; ik++) {
	    k = (ik-i0)*dk;
	    k *= k;
	    e = expf(-a*k);
	    e2 = e*e;
	    le = k*e;
	    
	    f = data[ik];
	    
	    ee += e2;
	    fe += f*e;
	    lfe += f*le;
	    le2 += le*e;
	    l2e2 += le*le;
	}
              
	aa = fe/ee;
	da = (2.*aa*le2 - lfe)/ee;
	num = aa*((aa*le2 - lfe) + eps*da);
	den = aa*aa*l2e2 + lfe/ee*(lfe - 2.*aa*le2) + eps*(1.-da*da);
        
	da = num/den;
        
	r2 = f2 - fe*fe/ee;
	if (verb) sf_warning("monof: iter=%d r2=%g da=%g aa=%g a=%g",
			     iter,r2,da,aa,a);

	a += da;
	if (r2 < eps || da*da < eps) break;
    }
        
    for (ik = 0; ik < nk; ik++) {
	k = (ik-i0)*dk;
	k *= k;
	data[ik] = aa*exp(-a*k);
    }
 
    if (verb) sf_warning ("%d iterations", iter);

    return a;
}

/* 	$Id: monof.c,v 1.1 2004/02/14 06:57:16 fomels Exp $	 */
