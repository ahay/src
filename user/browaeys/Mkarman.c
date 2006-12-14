/* Estimating the logarithm of a von Karman 1-D spectrum by nonlinear least squares */
/*
  Copyright (C) 2006 University of Texas at Austin
  
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
  Input data and function

  f  = log(data)
  l  = log(1 + a*k*k)
  aa = -(nu/2+1/4)

  Formulas for separable least squares:

  aa = f.l/(l.l + eps)
  ap = (f.lp-2*a*l.lp)/(l.l + eps)
  num = a*(a*l.lp + ap*(l.l + 2.*eps))
  den = ap*ap*l.l + 2*a*ap*l.lp + a*a*lp.lp
  da = num/den

  Requires:

  fl   -> f.l
  ll   -> l.l + eps
  flp  -> f.lp
  llp  -> l.lp
  lplp -> lp.lp
*/


int main(int argc, char* argv[])
{

    float *data /*input [nk] */; 
    int niter   /* number of iterations */; 
    float a0    /* initial value */;
    int nk      /* data length */;
    float dk    /* data sampling */; 
    float k0    /* initial location */;
    bool verb   /* verbosity flag */;

    int ik, iter;
    float f2, fl, ll, flp, llp, lplp, k, da, aa, f, l2;
    float lp, eps, num, den, r2, a, l, s;
    sf_file in, out;

    /* Estimate shape (Caution: data gets corrupted) */ 

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nk)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dk)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&k0)) sf_error("No o1= in input");
 
    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */
    if (!sf_getfloat("a0",&a0)) a0=-0.25;
    /* initial slope value */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    data = sf_floatalloc(nk);

    eps = 10.*FLT_EPSILON;
    eps *= eps;
    
    sf_floatread(data,nk,in);
    f2 = 0.;
    for (ik=0; ik < nk; ik++) {
	f = log(data[ik]);
	f2 += f*f;
    }

    a = a0; /* initial a */
    aa = 0.25;

    if (verb) sf_warning("got a0=%g k0=%g niter=%d nk=%d dk=%g",
			 a0,k0,niter,nk,dk);
	
    /* Gauss-Newton iterations */
    for (iter = 0; iter < niter; iter++) { 
	ll = eps;
	fl = flp = llp = lplp = 0.;
	for (ik = 0; ik < nk; ik++) {
	    k = ik*dk + k0;
	    k *= k;
            s = 1 + a*k;
	    l = log(s);
	    l2 = l*l;
	    lp = k/s; /* derivative of l with respect to a */
	    
	    f = log(data[ik]);
	    
	    ll += l2;
	    fl += f*l;
	    flp += f*lp;
	    llp += lp*l;
	    lplp += lp*lp;
	}
              
	aa = fl/ll;  /* amplitude */
	da = (flp - 2.*aa*llp)/ll;
	num = aa*(aa*llp + da*(ll+eps));
	den = aa*aa*lplp + da*(2.*aa*llp + da*(ll-eps));
        
	da = num/den; /* delta a */
        
	r2 = f2 - aa*aa*(ll+eps); /* residual squared */
	if (verb) sf_warning("iter=%d r2=%g da=%g aa=%g a=%g",
			     iter,r2,da,aa,a);

	a += da;     /* update a */
	if (r2 < eps || da*da < eps) break;
    }
        
    for (ik = 0; ik < nk; ik++) {
	k = ik*dk + k0;
	k *= k;
	data[ik] = exp(aa*log(1+a*k));
    }
 
    if (verb) sf_warning ("%d iterations", iter);        
	

    /* 
    Optimized parameters for f = log(data) = aa*log(1+a*k*k)
    with a=b*b and aa = -(nu/2+1/4)
    */
    sf_warning ("b=%g nu=%g",sqrt(a),-2*aa-0.5);
    sf_floatwrite (data,nk,out);
    
    exit (0);
}

/* 	$Id$	 */
