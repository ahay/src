#include <float.h>
#include <math.h>

#include <rsf.h>

#include "monof.h"

/* formulas for separable least squares:

a = f.e/(e.e + eps)
ap = (f.ep-2*a*e.ep)/(e.e + eps)
num = a*(a*e.ep + ap*(e.e + 2.*eps))
den = ap*ap*e.e + 2*a*ap*e.ep + a*a*ep.ep
da = num/den

need: 

fe -> f.e
ee -> e.e + eps
fep -> f.ep
eep -> e.ep
epep -> ep.ep

*/

/* destroys data */
float monof(float *data, int i0, int niter, float a0, int nk, 
	    float dk, bool verb)
{
    int ik, iter;
    float f2, fe, ee, fep, eep, epep, k, da, aa, f, e2;
    float ep, eps, num, den, r2, a, e;
    
    eps = 10.*FLT_EPSILON;
    eps *= eps;

    f2 = 0.;
    for (ik=0; ik < nk; ik++) {
	f = data[ik];
	f2 += f*f;
    }

    a = a0;
    aa = 1.;

    if (verb) sf_warning("got a0=%g i0=%d niter=%d nk=%d dk=%g",
			 a0,i0,niter,nk,dk);
	
    for (iter = 0; iter < niter; iter++) {
	ee = eps;
	fe = fep = eep = epep = 0.;
	for (ik = 0; ik < nk; ik++) {
	    k = (ik-i0)*dk;
	    k *= k;
	    e = expf(-a*k);
	    e2 = e*e;
	    ep = -k*e;
	    
	    f = data[ik];
	    
	    ee += e2;
	    fe += f*e;
	    fep += f*ep;
	    eep += ep*e;
	    epep += ep*ep;
	}
              
	aa = fe/ee;
	da = (fep - 2.*aa*eep)/ee;
	num = aa*(aa*eep + da*(ee+eps));
	den = aa*aa*epep + da*(2.*aa*eep + da*(ee-eps));
        
	da = num/den;
        
	r2 = f2 - aa*aa*(ee+eps); /* residual squared */
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

/* 	$Id: monof.c,v 1.2 2004/02/25 16:15:10 fomels Exp $	 */
