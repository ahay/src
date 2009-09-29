/* Estimation of von Karman autocorrelation 2D spectrum by nonlinear separable least squares. */
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
  Input in 2D:

  f  = log(data)
  l  = log(1 + k.A.k) where k = (kx ,ky) and A real symmetric correlation matrix
  aa = -(nu/2+1/2)
  lp = derivative of l with respect to A

  Formulas for nonlinear separable least squares and Gauss Newton:

  aa = f.l/(l.l + eps)                          linear slope 
  ap = (f.lp-2*aa*l.lp)/(l.l + eps)             derivative of aa with respect to A
  num = aa*(aa*l.lp + ap*(l.l + 2.*eps))
  den = ap*ap*l.l + 2*aa*ap*l.lp + aa*aa*lp.lp
  da = num/den                                  increment for A

  Requires:

  fl   -> f.l
  ll   -> l.l + eps
  flp  -> f.lp
  llp  -> l.lp
  lplp -> lp.lp
*/

int main(int argc, char* argv[])
{

    float **data               /* input [ny][nx] */;
    int niter                  /* number of iterations */; 
    int nx, ny                 /* axis */;
    float dx, dy, x0, y0       /* axis and initials */;
    bool verb                  /* verbosity flag */;

    int ix, iy, iter, j;
    float f2, fl, ll, flp[3], llp[3], lplp[6], x, y, da[3], aa, f, l2;
    float lp[3], eps, num[3], den[6], r2, a[3], l, x2, y2, xy, det, r;
    sf_file in, out;
    
    /* Estimate shape (Caution: data gets corrupted) */ 

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&ny)) sf_error("No n2= in input");

    if (!sf_histfloat(in,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&x0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d2",&dy)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&y0)) sf_error("No o2= in input");

    if (!sf_getfloat("a0",&a[0])) a[0]=1000.;
    /* starting correlation length in xx */
    if (!sf_getfloat("b0",&a[1])) a[1]=0.;
    /* starting correlation length in xy */
    if (!sf_getfloat("c0",&a[2])) a[2]=400.;
    /* starting correlation length in yy */

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    data = sf_floatalloc2(nx,ny);

    /* Inversion */

    eps = 10.*FLT_EPSILON;
    eps *= eps;
	    
    sf_floatread(data[0],nx*ny,in);

    f2 = 0.;
    for (iy=0; iy < ny; iy++) {
        for (ix=0; ix < nx; ix++) {
            f = log(data[iy][ix]);
	    f2 += f*f;		    
	}
    }
    aa = -0.5;

    if (verb) sf_warning("got a0=%g b0=%g c0=%g niter=%d\n"
		         "nx=%d dx=%g x0=%g ny=%d dy=%g y0=%g",
		         a[0],a[1],a[2],niter,nx,dx,x0,ny,dy,y0);

    /* Gauss-Newton iterations */

    for (iter = 0; iter < niter; iter++) {
        ll = eps;
        fl = 0.;
        for (j=0; j < 3; j++) {	    
         	flp[j] = 0.;
		llp[j] = 0.;
		lplp[j] = 0.;
	        lplp[3+j] = 0.;
        }
        for (iy=0; iy < ny; iy++) {
        	y = y0+iy*dy;
		y2 = y*y;
		for (ix=0; ix < nx; ix++) {
		    x = x0+ix*dx;
		    x2 = x*x;
		    xy = x*y;
                    r = 1 + a[0]*x2 + a[1]*xy + a[2]*y2;
	            l = log(r);
	  	    l2 = l*l;

	            /* Derivative of l with respect to a */

		    lp[0] = x2/r;
		    lp[1] = xy/r;
		    lp[2] = y2/r;
	    
		    f = log(data[iy][ix]);
	    
		    ll += l2;
		    fl += f*l;
		    for (j=0; j < 3; j++) {
			flp[j] += f*lp[j];
			llp[j] += l*lp[j];
		    }
		    lplp[0] += lp[0]*lp[0];
		    lplp[1] += lp[0]*lp[1];
		    lplp[2] += lp[0]*lp[2];
		    lplp[3] += lp[1]*lp[1];
		    lplp[4] += lp[1]*lp[2];
		    lplp[5] += lp[2]*lp[2];
		}
	}
           
     	aa = fl/ll;  /* amplitude */

	for (j=0; j < 3; j++) {
		da[j] = (flp[j] - 2.*aa*llp[j])/ll;
		num[j] = aa*(aa*llp[j] + da[j]*(ll+eps));
        }
        den[0] = aa*aa*lplp[0] + da[0]*(2.*aa*llp[0] + da[0]*(ll-eps));
        den[1] = aa*(aa*lplp[1] + da[1]*llp[0]) + 
	       da[0]*(aa*llp[1] + da[1]*(ll-eps));
	den[2] = aa*(aa*lplp[2] + da[2]*llp[0]) + 
	       da[0]*(aa*llp[2] + da[2]*(ll-eps));
	den[3] = aa*aa*lplp[3] + da[1]*(2.*aa*llp[1] + da[1]*(ll-eps));
        den[4] = aa*(aa*lplp[4] + da[1]*llp[2]) + 
	       da[2]*(aa*llp[1] + da[1]*(ll-eps));
	den[5] = aa*aa*lplp[5] + da[2]*(2.*aa*llp[2] + da[2]*(ll-eps));
	    
        det = 
	    den[2]*(den[2]*den[3] - den[1]*den[4]) + 
	    den[1]*(den[1]*den[5] - den[2]*den[4]) + 
	    den[0]*(den[4]*den[4] - den[3]*den[5]);
	    
        if (det > 0. && eps > det) {
		det = eps;
        } else if (det < 0. && -eps < det) {
		det = -eps;
        }
	
	da[0] = (
		num[2]*(den[2]*den[3] - den[1]*den[4]) + 
		num[1]*(den[1]*den[5] - den[2]*den[4]) + 
		num[0]*(den[4]*den[4] - den[3]*den[5])
		)/det;
	da[1] = ( 
		num[1]*(den[2]*den[2] - den[0]*den[5]) + 
		num[0]*(den[1]*den[5] - den[2]*den[4]) + 
		num[2]*(den[4]*den[0] - den[1]*den[2])
		)/det;
        da[2] = (
		num[0]*(den[2]*den[3] - den[1]*den[4]) + 
		num[2]*(den[1]*den[1] - den[3]*den[0]) + 
		num[1]*(den[4]*den[0] - den[1]*den[2])
		)/det;
	    
	r2 = f2 - aa*aa*(ll+eps); /* residual squared */

        if (verb) sf_warning("iter=%d r2=%g da=(%g,%g,%g) "
				 "aa=%g a=(%g,%g,%g)",
				 iter,r2,da[0],da[1],da[2],aa,a[0],a[1],a[2]);
  	
        /* Update a */
	for (j=0; j < 3; j++) {
		a[j] += da[j];
        }
	if (a[0] < 0.) a[0] = 0.;
        if (a[2] < 0.) a[2] = 0.;
	    
        if (r2 < eps || da[0]*da[0] + da[1]*da[1] + da[2]*da[2] < eps) break;
    } 
    /* iter */
    
    for (iy=0; iy < ny; iy++) {
	    y = y0+iy*dy;
	    y2 = y*y;
	    for (ix=0; ix < nx; ix++) {
		x = x0+ix*dx;
		x2 = x*x;
		xy = x*y;
		data[iy][ix] =  exp(aa*log(1+a[0]*x2+a[1]*xy+a[2]*y2));
	    }
    }


    if (verb) sf_warning ("%d iterations", iter);

    /* 
         Optimized parameters for f = log(data) = aa*log(1+a[0]*x2+a[1]*xy+a[2]*y2)
         with a=b*b and aa = -(nu/2+1/2)
    */

    sf_warning ("axx=%g axy=%g ayy=%g nu=%g",a[0],a[1],a[2],-2*aa-1.);
    
    sf_floatwrite (data[0],nx*ny,out);

    exit (0);
}

/* 	$Id$	 */
