/* Inversion for von Karman autocorrelation 1D spectrum. */
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

/*
  Input data and functions:

  f  = log[data(k)]
  l  = log[1+b*b*k*k]
  aa = -(nu/2+1/4)
  dd = log(d0)

  s  = data length
  sf = sum of f  
  sl = sum of l

  lp = derivative of l with respect to b
  sp = derivative of sl with respect to b

  Formulae for nonlinear separable least squares for parameters (aa,dd)

  aa = [s*f.l-(sl+eps)*sf]/[s*l.l-(sl+eps)^2]      linear slope 
  dd = [sf*l.l-(sl+eps)*f.l]/[s*l.l-(sl+eps)^2]    origin constant

  Derivatives with respect to b
 
  ap = [s*f.lp-sp*sf-2*aa*(s*l.lp-(sl+eps)*sp)]/[s*l.l-(sl+eps)^2]                   linear slope
  dp = [2*sf*l.lp-sp*f.l-(sl+eps)*f.lp-2*dd*(s*l.lp-(sl+eps)*sp)]/[s*l.l-(sl+eps)^2] origin constant

  Gauss Newton inversion for nonlinear parameter b

  num = ap*f.l + aa*f.lp + dp*sf
      - aa*(ap*l.l + aa*l.lp + dp*(sl+eps))
      - dd*(ap*(sl+eps) + aa*sp + dp*s)   
  den = ap*ap*l.l + aa*aa*lp.lp + dp*dp*s
      + 2*(aa*ap*l.lp + dp*ap*(sl+eps) + dp*aa*sp)
  db = num/den

  Requires:

  s,sf,lp,sp

  sl   -> sl + eps
  fl   -> f.l
  ll   -> l.l
  flp  -> f.lp
  llp  -> l.lp
  lplp -> lp.lp
*/


int main(int argc, char* argv[])
{

    float *data /*input [nk] */; 
    int niter   /* number of iterations */; 
    float b0    /* initial value */;
    int nk      /* data length */;
    float dk    /* data sampling */; 
    float k0    /* initial location */;
    bool verb   /* verbosity flag */;

    int ik, iter;
    float f, l, fl, ll, flp, llp, lplp, s, sf, sl, lp, sp;
    float aa, dd, ap, dp, num, den, eps;
    float k, r, f2, l2, b, db, r2, vv, dv;

    sf_file in, out, pa;

    /* Estimate shape (Caution: data gets corrupted) */ 

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    pa = sf_output("pa");


    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nk)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dk)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&k0)) sf_error("No o1= in input");
 
    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */
    if (!sf_getfloat("b0",&b0)) b0=50.;
    /* initial nonlinear parameter value */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    sf_putint(pa,"n1",3);
    sf_putint(pa,"nk",nk);
    sf_putfloat(pa,"dk",dk);
    sf_putfloat(pa,"k0",k0);
    sf_fileflush(pa,in);

    data = sf_floatalloc(nk);

    eps = 10.*FLT_EPSILON;
    eps *= eps;
    
    sf_floatread(data,nk,in);
    f2 = sf = s = 0.;
    for (ik=0; ik < nk; ik++) {
	f = log(data[ik]);
	f2 += f*f;
        sf += f;
        s += 1.;
    }

    b = b0; /* initial b */
    aa = -1.;
    dd = 1.;

    if (verb) sf_warning("got b0=%g k0=%g niter=%d nk=%d dk=%g",
			 b0,k0,niter,nk,dk);
	
    /* Gauss-Newton iterations */
    for (iter = 0; iter < niter; iter++) {
        sl = eps;
	ll = fl = flp = llp = lplp = sp = 0.;
	for (ik = 0; ik < nk; ik++) {
	    k = ik*dk + k0;
	    k *= k;
            r = 1 + b*b*k;
	    l = log(r);
	    l2 = l*l;
	    lp = 2.*b*k/r; /* derivative of l with respect to b */
	    sl += l;
	    sp += lp;

	    f = log(data[ik]);            
	    
	    ll += l2;
	    fl += f*l;
	    flp += f*lp;
	    llp += lp*l;
	    lplp += lp*lp;
	}

        vv = s*ll - sl*sl;
	aa = (s*fl - sf*sl)/vv;  /* amplitude slope */
	dd = (sf*ll - sl*fl)/vv; /* amplitude constant */

        dv = s*llp - sl*sp;
        ap = (s*flp - sf*sp - 2.*aa*dv)/vv;              /* derivative slope */
        dp = (2.*sf*llp - sp*fl - sl*flp - 2.*dd*dv)/vv; /* derivative constant */

        num = ap*fl + aa*flp + dp*sf - aa*(ap*ll + aa*llp + dp*sl) - dd*(sl*ap + sp*aa + s*dp);
        den = ap*ap*ll + aa*aa*lplp + s*dp*dp + 2.*(aa*ap*llp + sl*dp*ap + sp*dp*aa); 

	db = num/den; /* delta b */

	r2 = f2 - (s*dd*dd + aa*aa*ll + 2.*dd*aa*(sl+eps)); /* residual squared */
	if (verb) sf_warning("iter=%d r2=%g db=%g aa=%g dd=%g b=%g",
			     iter,r2,db,aa,dd,b);

	b += db;     /* update b */
	if (r2 < eps || db*db < eps) break;
    }
        
    for (ik = 0; ik < nk; ik++) {
	k = ik*dk + k0;
	k *= k;
	data[ik] = exp(dd+aa*log(1+b*b*k));
    }
 
    if (verb) sf_warning ("%d iterations", iter);        
	
    /* Optimized parameters for f = log(data) = dd + aa*log(1+b*b*k*k) where aa = -(nu/2+1/4) */
    sf_warning ("b=%g nu=%g d0=%g",b,-2*aa-0.5,exp(dd));

    sf_floatwrite(&b,1,pa);
    sf_floatwrite(&aa,1,pa);
    sf_floatwrite(&dd,1,pa);

    sf_floatwrite (data,nk,out);
    
    exit (0);
}

/* 	$Id$	 */
