/* Ricker wavelet estimation. */
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

int main(int argc, char* argv[])
{
    int n2, i2, na, ia, niter, iter;
    float f0, df, f, f2, m, m0, m2, m3, *data, *r, *rp;
    float rd, r2, rpd, rp2, rpr, ap, num, den, a, e, dm, eps, di, d2;
    bool verb;
    sf_file in, out, ma;

    sf_init (argc,argv);
    in  = sf_input("in");
    out = sf_output("out");
    ma  = sf_output("ma");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&na)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_histfloat(in,"d1",&df)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&f0)) sf_error("No o1= in input");

    if (!sf_getfloat("m",&m0)) m0=f0+0.25*(na-1)*df;
    /* initial frequency */
    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    sf_putint(ma,"n1",2);
    sf_putint(ma,"nf",na);
    sf_putfloat(ma,"df",df);
    sf_putfloat(ma,"f0",f0);
    sf_fileflush(ma,in);

    data = sf_floatalloc(na);
    r = sf_floatalloc(na);
    rp = sf_floatalloc(na);

    eps = 10.*FLT_EPSILON;
    eps *= eps;

    for (i2=0; i2 < n2; i2++) {
	sf_warning("slice %d of %d;",i2+1,n2);
	sf_floatread(data,na,in);

	d2 = 0.;
	for (ia=0; ia < na; ia++) {
	    di = data[ia];
	    if (di < 0.) data[ia] = -di;
	    d2 += di*di;
	}

	m = m0;

	for (iter = 0; iter < niter; iter++) {
	    m2 = m*m;
	    m3 = m2*m;
           
	    rd = r2 = rpd = rp2 = rpr = 0.;
	    for (ia = 0; ia < na; ia++) {
		f = f0 + ia*df;
		f2 = f*f;
		e = exp(-f2/m2);
              
		r[ia] = e*f2/m2;
		rp[ia] = 2.*e*f2*(f2-m2)/(m3*m2); /* dr/dm */
		
		rd += r[ia]*data[ia];
		r2 += r[ia]*r[ia];
		rpd += rp[ia]*data[ia];
		rp2 += rp[ia]*rp[ia];
		rpr += rp[ia]*r[ia];
	    }
              
	    a = rd/(r2 + eps);
	    ap = (rpd-2.*rpr*a)/(r2 + eps); /* da/dm */

	    num =  a*(rpd-rpr*a)+ap*(rd-r2*a);
	    den = a*a*rp2 + 2.*a*ap*rpr + ap*ap*r2 + eps;
        
	    dm = num/den;
        
	    r2 = d2 - 2.*rd*a + r2*a*a; /* ||d - a*r||^2 */ 
	    rp2 = dm*dm;

	    if (verb && 5000 > n2) sf_warning("iter=%d r2=%g rp2=%g m=%g a=%g",
					      iter,r2,rp2,m,a);

	    m += dm;
	    if (r2 < eps || rp2 < eps) break;
	}
        
	m = fabsf(m);
	m2 = m*m;
        
	sf_floatwrite(&m2,1,ma);
	sf_floatwrite(&a,1,ma);
        
	for (ia = 0; ia < na; ia++) {
	    f = f0 + ia*df;
	    f2 = f*f;
	    data[ia] = a*exp(-f2/m2)*f2/m2;
	}
        
	if (verb) sf_warning("m=%g a=%g",m,a*m*sqrtf(SF_PI)*0.5);
	if (verb) sf_warning ("%d of %d, %d iterations", i2+1, n2, iter);
        
	sf_floatwrite (data,na,out);
    }

    sf_warning(".");

    exit (0);
}

/* 	$Id$	 */
