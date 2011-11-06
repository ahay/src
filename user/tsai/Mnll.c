/* Nonlinear optimization for von Karman autocorrelation 1D spectrum with weighted least squares */
/* 1 - Separable least squares for 2 parameters */
/* 2 - Gauss Newton on 1 nonlinear parameter */
/*
  Copyright (C) 2007 University of Texas at Austin
  
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

#include "karman.h"

/*
  Input model

  f = log[data(k)]
  l = log[1+x[0]*k*k]

  x[0] = squared correlation length
  y = -(Hu/2+1/4)
  z = log(amplitude)

  s  = data length
  sf = sum of f  
  sl = sum of l

  lp = derivative of l with respect to x[0]
  sp = derivative of sl with respect to x[0]

  Formulae for nonlinear separable least squares for parameters (y,z)

  y = [s*f.l-(sl+eps)*sf]/[s*l.l-(sl+eps)^2]      linear slope y
  z = [sf*l.l-(sl+eps)*f.l]/[s*l.l-(sl+eps)^2]    origin constant z

  Derivatives with respect to x[0]
 
  yp = [s*f.lp-sp*sf-2*y*(s*l.lp-(sl+eps)*sp)]/[s*l.l-(sl+eps)^2]
  zp = [2*sf*l.lp-sp*f.l-(sl+eps)*f.lp-2*z*(s*l.lp-(sl+eps)*sp)]/[s*l.l-(sl+eps)^2]

  Gauss Newton iteration for nonlinear parameter x[0]

  num = yp*f.l + y*f.lp + zp*sf
      - y*(yp*l.l + y*l.lp + zp*(sl+eps))
      - z*(yp*(sl+eps) + y*sp + zp*s)   
  den = yp*yp*l.l + y*y*lp.lp + zp*zp*s
      + 2*(y*yp*l.lp + zp*yp*(sl+eps) + zp*y*sp)
  dx = num/den

  Requires:

  sf,lp,sp

  s    -> s + eps
  sl   -> sl + eps
  ll   -> l.l + eps
  fl   -> f.l
  ll   -> l.l
  flp  -> f.lp
  llp  -> l.lp
  lplp -> lp.lp
*/

void karman(float *data                /* input [nk] */,
	    float *pred                /* prediction [nk] */,
	    int nliter                 /* number of reweighting iterations */,
	    int niter                  /* number of iterations */, 
	    float* x                   /* estimated parameters [3] */, 
	    int nk, float dk, float k0 /* axis */, 
	    bool verb                  /* verbosity flag */)

            /*< Caution: data gets corrupted >*/ 

{
    int ik, iter, i;
    float f, l, fl, ll, flp, llp, lplp, s, sf, sl, lp, sp;
    float y, z, yp, zp, num, den, eps, dx, k, r, f2, w, vv, dv, r2;

    eps = 10.*FLT_EPSILON;
    eps *= eps;

    dk = 2.*SF_PI*dk;
    k0 = 2.*SF_PI*k0;

    if (verb) sf_warning("got x0=%g y0=%g z0=%g\n"
			 "nk=%d dk=%g k0=%g",
			 x[0],x[1],x[2],nk,dk,k0);

    if (nliter > 1) sf_irls_init(nk);

    sf = 0.;
    for (ik=0; ik < nk; ik++) {
        sf += log(data[ik]);
    }

    x[0] *= x[0];
    y = x[1];
    z = x[2];

    for (i = 0; i < nliter; i++) {

	f2 = 0.;
	for (ik=0; ik < nk; ik++) {
	    f = log(data[ik]);
	    w = (0 == i)? 1.:pred[ik];
	    f2 += w*f*f;		    
	}

	/* Gauss Newton iterations */
	s = (float)nk + eps;
        for (iter = 0; iter < niter; iter++) {
	    sl = ll = eps;
	    fl = flp = llp = lplp = sp = 0.;
	    for (ik = 0; ik < nk; ik++) {
		k = ik*dk + k0;
		k *= k;
		r = 1. + x[0]*k;
		l = log(r);
	        lp = k/r /* derivative of l with respect to x */;

		f = log(data[ik]);            
	    	w = (0 == i)? 1.:pred[ik];
				
		sl += w*l;
		sp += w*lp;
		ll += w*l*l;
		fl += w*f*l;
		flp += w*f*lp;
		llp += w*lp*l;
		lplp += w*lp*lp;
	    }

            vv = 1./(s*ll - sl*sl);
	    y = (s*fl - sf*sl)*vv; 
	    z = (sf*ll - sl*fl)*vv; 

	    dv = s*llp - sl*sp;
	    yp = (s*flp - sf*sp - 2.*y*dv)*vv;
	    zp = (2.*sf*llp - sp*fl - sl*flp - 2.*z*dv)*vv;

	    num = yp*fl + y*flp + zp*sf - y*(yp*ll + y*llp + zp*sl) - z*(sl*yp + sp*y + s*zp);
	    den = yp*yp*ll + y*y*lplp + s*zp*zp + 2.*(y*yp*llp + sl*zp*yp + sp*zp*y);

	    dx = num/den /* delta x */;

	    r2 = f2 - ((s+eps)*z*z + y*y*(ll+eps) + 2.*z*y*(sl+eps)) /* residual squared */;
	    if (verb) sf_warning("iter=%d r2=%g dx=%g x=%g y=%g z=%g",
			          iter,r2,dx,sqrt(x[0]),-2.*y-0.5,exp(z));
	    x[0] += dx /* update x */;
	    if (r2 < eps || dx*dx < eps) break;
        } /* iter */

        for (ik = 0; ik < nk; ik++) {
	    k = ik*dk + k0;
	    k *= k;
	    l = z + y*log(1.+x[0]*k);
	    pred[ik] = (i < nliter-1)? l-data[ik]: exp(l);
	}

	if (i < nliter-1) sf_cauchy (nk,pred,pred);
	
	if (verb) sf_warning ("%d iterations", iter);
    } /* reweighting iterations */

    if (nliter > 1) sf_irls_close();
    /* Optimized parameters for f = log(data) = z + y*ln(1+x*k*k) */
    sf_warning ("x=%g Hu=%g A0=%g",sqrt(x[0]),-2.*y-0.5,exp(z));
}

/* 	$Id$	 */

