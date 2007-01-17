/* Nonlinear optimization for von Karman autocorrelation 1D spectrum. */
/* 1. Separable least squares for 2 parameters */
/* 2. Full Newton on 1 non linear parameter */
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
  Input model

  f = ln[data(k)]
  x = squared correlation length

  s = 1 + x*k*k
  l = ln(s)
   
  r(k,x) = z(x) + y(x)*l(k,x) - f(k) = residual component

  Derivatives of l with respect to x

  lp = k*k/s
  ls = -lp*lp

  Separable least squares

  y = -(Hu/2 + 1/4)
  z = ln[amplitude]
  
  y = (nk*fl - sl*sf)*dv
  z = (sf*ll - sl*fl)*dv

  yp = (nk*flp - slp*sf - y*dp)*dv
  zp = (2.*sf*llp - slp*fl - sl*flp - z*dp)*dv

  ys = (nk*fls - sls*sf - y*ds - 2.*yp*dp)*dv
  zs = (2.*(sf*(lls-sls) - slp*flp - zp*dp) - sls*fl - sl*fls - z*ds)*dv
  
  where

  dv = 1./(nk*ll - sl*sl)
  dp = 2.*(nk*llp - sl*slp)
  ds = 2.*(nk*(lls-sls) - sl*sls - slp*slp)

  Full Newton method on non linear parameter x

  num = zp*(nk*z + y*sl - sf) + yp*(z*sl + y*ll - fl) + y*(z*slp + y*llp - flp)

  sj = zp*zp*nk + yp*yp*ll - y*sls + 2.*[zp*(yp*sl+y*slp) + yp*y*llp]
  sp = zs*(nk*z - sf) + sl*(z*ys + y*zs) + 2.*yp*(z*slp + y*llp - flp) + ys*(y*ll - fl) + y*(y*lls - z*sls - fls)

  dx = - num/(sj+sp)

  Requires

  sl = sum of l
  sf = sum of f

  slp = sum of lp
  sls = sum of ls

  ll = sum of l*l
  fl = sum of f*l

  llp = sum of l*lp
  lls = sum of l*ls

  flp = sum of f*lp
  fls = sum of f*ls
*/


int main(int argc, char* argv[])
{

    float *data /*input [nk] */; 
    int niter   /* number of iterations */; 
    float x0    /* initial value */;
    int nk      /* data length */;
    float dk    /* data sampling */; 
    float k0    /* initial location */;
    bool verb   /* verbosity flag */;

    int ik, iter;
    float k, f, x, s, l, lp, ls;
    float m, sl, sf, ll, fl, slp, llp, flp, sls, fls, lls;
    float dx, sj, sp, num ,dv, dp, ds, y, z, yp, zp, ys, zs, r, r2, eps;
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
    if (!sf_getfloat("x0",&x0)) x0=1.;
    /* initial non linear parameter value */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    data = sf_floatalloc(nk);
    sf_floatread(data,nk,in);

    x = x0; /* initial x */

    eps = 10.*FLT_EPSILON;
    eps *= eps;

    if (verb) sf_warning("got x0=%g k0=%g niter=%d nk=%d dk=%g",
			 x0,k0,niter,nk,dk);
	
    /* Full Newton iterations */
    m = (float)nk;
    for (iter = 0; iter < niter; iter++) {
        sl = sf = ll = fl = slp = llp = flp = sls = fls = lls = r2 = 0.;
	for (ik = 0; ik < nk; ik++) {
	    k = ik*dk + k0;
	    k *= k;
            s = 1. + x*k;
	    l = log(s);
	    lp = k/s; /* derivative of l with respect to x */
            ls = -lp*lp;

	    sl += l;
	    ll += l*l;
	    slp += lp;
	    llp += l*lp;
	    sls += ls;
	    lls += l*ls;
	    
	    f = log(data[ik]);
	    r = z + y*l - f;

            r2 += r*r;

	    sf += f;
	    fl += f*l;
	    flp += f*lp;
	    fls += f*ls;
	}

	dv = 1./(m*ll - sl*sl);
	y = (m*fl - sl*sf)*dv;
	z = (sf*ll - sl*fl)*dv;

	dp = 2.*(m*llp - sl*slp);
	yp = (m*flp - slp*sf - y*dp)*dv;
	zp = (2.*sf*llp - slp*fl - sl*flp - z*dp)*dv;
        ds = 2.*(m*(lls-sls) - sl*sls - slp*slp);
	ys = (m*fls - sls*sf - 2.*yp*dp - y*ds)*dv;
	zs = (2.*(sf*(lls-sls) - slp*flp - zp*dp) - sls*fl - sl*fls - z*ds)*dv;
  
	num = zp*(m*z + y*sl - sf) + yp*(z*sl + y*ll - fl) + y*(z*slp + y*llp - flp);
	sj = zp*zp*m + yp*yp*ll - y*sls + 2.*(zp*(yp*sl+y*slp) + yp*y*llp);
	sp = zs*(m*z - sf) + sl*(z*ys + y*zs) + 2.*yp*(z*slp + y*llp - flp) 
           + ys*(y*ll - fl) + y*(y*lls - z*sls - fls);

	dx = - num/(sj+sp);

        /* Residual squared */
	for (ik = 0; ik < nk; ik++) {
	    k = ik*dk + k0;
	    k *= k;
	    f = log(data[ik]);
	    r = z + y*log(1.+x*k) - f;

            r2 += r*r;
	}
	
	if (verb) sf_warning("iter=%d r2=%g dx=%g x=%g y=%g z=%g",
			     iter,r2,dx,x,y,z);

	if (r2 < eps || dx*dx < eps) break;

	x += dx; /* update x */
    }
        
    for (ik = 0; ik < nk; ik++) {
	k = ik*dk + k0;
	k *= k;
	data[ik] = exp(z + y*log(1.+x*k));
    }
 
    if (verb) sf_warning ("%d iterations", iter);        
	

    /* Optimized parameters for f = log(data) = z + y*ln(1+x*k*k) */
    sf_warning ("x=%g nu=%g",sqrt(x),-2.*y-0.5);
    sf_floatwrite (data,nk,out);
    
    exit (0);
}

/* 	$Id$	 */
