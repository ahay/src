/* Fitting a Gaussian to 1-D spectrum by nonlinear least squares */
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

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

float monof(float *data /*input [nk] */, 
	    int i0      /* maximum location */, 
	    int niter   /* number of iterations */, 
	    float a0    /* initial value */, 
	    int nk      /* data length */, 
	    float dk    /* data sampling */, 
	    bool verb   /* verbosity flag */)
/*< Estimate shape. (Caution: data gets corrupted) >*/ 
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

    a = a0; /* initial a */
    aa = 1.;

    if (verb) sf_warning("got a0=%g i0=%d niter=%d nk=%d dk=%g",
			 a0,i0,niter,nk,dk);
	
    /* Gauss-Newton iterations */
    for (iter = 0; iter < niter; iter++) { 
	ee = eps;
	fe = fep = eep = epep = 0.;
	for (ik = 0; ik < nk; ik++) {
	    k = (ik-i0)*dk;
	    k *= k;
	    e = expf(-a*k);
	    e2 = e*e;
	    ep = -k*e; /* derivative of e with respect to a */
	    
	    f = data[ik];
	    
	    ee += e2;
	    fe += f*e;
	    fep += f*ep;
	    eep += ep*e;
	    epep += ep*ep;
	}
              
	aa = fe/ee;  /* amplitude */
	da = (fep - 2.*aa*eep)/ee;
	num = aa*(aa*eep + da*(ee+eps));
	den = aa*aa*epep + da*(2.*aa*eep + da*(ee-eps));
        
	da = num/den; /* delta a */
        
	r2 = f2 - aa*aa*(ee+eps); /* residual squared */
	if (verb) sf_warning("monof: iter=%d r2=%g da=%g aa=%g a=%g",
			     iter,r2,da,aa,a);

	a += da;     /* update a */
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

/* 	$Id: monof.c 7107 2011-04-10 02:04:14Z ivlad $	 */
