/* Inversion for von Karman autocorrelation 1D spectrum. */
/* Full Newton on 3 parameters */

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
  Input data and model vector

  f  = ln[data(k)]

  x1 = ln(amplitude)
  x2 = -(nu/2+1/4)
  x3 = correlation length

  s = 1+x3*x3*k*k
  l = ln(s)

  Full Newton method with component of R(x) residual vector
   
  r(k) = x1 + x2*l(k,x3) - f

  Jacobian 

  J(k,1) = 1
  J(k,2) = ln(s) = l
  J(k,3) = 2*x2*x3*k*k/s

  Requires

  sm = sum of l
  s2 = sum of l*l
  y1 = sum of residual r
  y2 = sum of r*l
  y3 = sum of 2*x2*x3*r*k*k/s
  sk = sum of k*k/s
  s1 = sum of k*k*(x1+2*x2*l-f)/s
  s3 = sum of k*k*(2*x2*x3*x3*k*k+(1-x3*x3*k*k)*r)/(s*s)

  [J(x)^t R(x)]^t = (y1,y2,y3)

  a = J^t(x)J(x) + S(x) is

  nk    sm    2*x2*x3*sk
        s2    2*x3*s1
              2*x2*s3

  Determinant should not be zero
  d = a11*(a22*a33 - a23*a23) + 2*a12*a13*a23 - a33*a12*a12 - a22*a13*a13
  b = symmetric inverse matrix (to divide by d) is

  a22*a33-a23*a23     a13*a23-a12*a33     a12*a23-a13*a22
                      a11*a33-a13*a13     a12*a13-a23*a11
                                          a11*a22-a12*a12
  Update

  dx1 = -1/d*[b11*y1 + b12*y2 + b13*y3]
  dx2 = -1/d*[b12*y1 + b22*y2 + b23*y3]
  dx3 = -1/d*[b13*y1 + b23*y2 + b33*y3]
*/


int main(int argc, char* argv[])
{

    float *data /*input [nk] */; 
    int nk      /* data length */;
    float dk    /* data sampling */; 
    float k0    /* initial location */;
    bool verb   /* verbosity flag */;
    int niter   /* number of iterations */; 

    int ik, iter;
    float x1, x2, x3, dx1, dx2, dx3, x10, x20, x30, x4, x5;
    float k, f, s, l, r, p, q, s2, y1, y2, y3, sk, s1, s3;
    float a11, a12, a13, a22, a23, a33, a2, a3, a4, a5, a6, d;
    float b11, b12, b13, b22, b23, b33, eps, r2; 
    
    sf_file in, out, prm;

    /* Estimate shape (Caution: data gets corrupted) */ 

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    prm = sf_output("prm");


    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nk)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dk)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&k0)) sf_error("No o1= in input");
 
    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */
    if (!sf_getfloat("x10",&x10)) x10=6.;
    /* initial nonlinear parameter x1 value */
    if (!sf_getfloat("x20",&x20)) x20=-0.5;
    /* initial nonlinear parameter x2 value */
    if (!sf_getfloat("x30",&x30)) x30=200.;
    /* initial nonlinear parameter x3 value */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    sf_putint(prm,"n1",3);
    sf_putint(prm,"nk",nk);
    sf_putfloat(prm,"dk",dk);
    sf_putfloat(prm,"k0",k0);
    sf_fileflush(prm,in);

    data = sf_floatalloc(nk);
    
    sf_floatread(data,nk,in);

    x1 = x10; /* initial x1 */
    x2 = x20; /* initial x2 */
    x3 = x30; /* initial x3 */
    x4 = x3*x3;

    eps = 10.*FLT_EPSILON;
    eps *= eps;
    a11 = (float)nk;

    if (verb) sf_warning("got x10=%g x20=%g x30=%g k0=%g niter=%d nk=%d dk=%g",
			 x10,x20,x30,k0,niter,nk,dk);
	
    /* Gauss iterations */
    for (iter = 0; iter < niter; iter++) {
        x5 = 2.*x2*x3;
        s2 = y1 = y2 = y3 = sk = s1 = s3 = a12 = r2 = 0.;
	for (ik = 0; ik < nk; ik++) {
	    k = ik*dk + k0;
	    k *= k;
            s = 1. + x4*k;
	    l = log(s);
            q = 1./s;
	    p = k*q;
            
            s2 += l*l;
            sk += p;
            a12 += l;

	    f = log(data[ik]);       
	    r = x1 + x2*l - f;

	    y1 += r;
	    y2 += r*l;
	    y3 += r*p;
	    s1 += p*(r + x2*l);
            s3 += p*q*(2.*x2*x4*k + r*(1.-x4*k));

	    r2 += r*r;
	}

        y3 *= x5;

	a13 = x5*sk;
	a22 = s2;
        a23 = 2.*x3*s1;
	a33 = 2.*x2*s3;

	a2 = a12*a12;
	a3 = a12*a13;
	a4 = a13*a13;
        a5 = a23*a23;
        a6 = a22*a33;

        d = 1./(a11*(a6-a5) + 2.*a3*a23 - a33*a2 - a22*a4);

	b11 = a6-a5;
        b12 = a13*a23-a12*a33;
	b13 = a12*a23-a13*a22;
	b22 = a11*a33-a4;
	b23 = a3-a23*a11;
	b33 = a11*a22-a2;

	dx1 = -d*(b11*y1 + b12*y2 + b13*y3);
	dx2 = -d*(b12*y1 + b22*y2 + b23*y3);
	dx3 = -d*(b13*y1 + b23*y2 + b33*y3);

	x1 += dx1;
	x2 += dx2;
	x3 += dx3;
        x4 = x3*x3;
	
        if (verb) sf_warning("iter=%d r2=%g dx1=%g dx2=%g dx3=%g x1=%g x2=%g x3=%g",
			     iter,r2,dx1,dx2,dx3,x1,x2,x3);
	if (r2 < eps || (dx1*dx1+dx2*dx2+dx3*dx3) < eps) break;    
    }
 
    for (ik = 0; ik < nk; ik++) {
	k = ik*dk + k0;
	k *= k;
	data[ik] = exp(x1+x2*log(1.+x4*k));
    }
 
    if (verb) sf_warning ("%d iterations", iter);        
	
    /* Optimized parameters for f = log(data) = log(d0) + a*log(1+b*b*k*k) where a = -(nu/2+1/4) */
    sf_warning ("b=%g nu=%g d0=%g",x3,-2*x2-0.5,exp(x1));

    sf_floatwrite(&x3,1,prm);
    sf_floatwrite(&x2,1,prm);
    sf_floatwrite(&x1,1,prm);

    sf_floatwrite (data,nk,out);

    exit (0);
}

/* 	$Id$	 */
