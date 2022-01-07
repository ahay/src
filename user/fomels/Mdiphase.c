/* Derivative of local frequency. */
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
#include <math.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    int nh, n1,n2, i1,i2, i, n12, niter, dim, n[SF_MAX_DIM], rect[SF_MAX_DIM];
    float *trace, *hilb, *dtrace, *dhilb, *num, *den, *phase, a,b,c, mean, d1;
    float *traced, *hilbd, *dtraced, *dhilbd, *dnum, *dden, da, db;
    char key[6];
    bool hertz, verb;
    sf_file in, out, df;
	
    sf_init (argc,argv);
    df = sf_input("in");
    in = sf_input("sig");
    out = sf_output("out");
	
    if (SF_FLOAT   != sf_gettype(in)) sf_error("Need float input");
    
    dim = sf_filedims (in,n);
    n1 = n[0];
    n12 = 1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
	n12 *= n[i];
    }
    n2 = n12/n1;

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity */
	
    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
	
    trace = sf_floatalloc(n1);
    hilb = sf_floatalloc(n1);

    dtrace = sf_floatalloc(n1);
    dhilb = sf_floatalloc(n1);

    traced = sf_floatalloc(n1);
    hilbd = sf_floatalloc(n1);

    dtraced = sf_floatalloc(n1);
    dhilbd = sf_floatalloc(n1);
	
    num = sf_floatalloc(n12);
    den = sf_floatalloc(n12);
    phase = sf_floatalloc(n12);

    dnum = sf_floatalloc(n12);
    dden = sf_floatalloc(n12);
	
    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */
	
    if (!sf_getint("order",&nh)) nh=100;
    /* Hilbert transformer order */
    if (!sf_getfloat("ref",&c)) c=1.;
    /* Hilbert transformer reference (0.5 < ref <= 1) */
	
    if (!sf_getbool("hertz",&hertz)) hertz=false;
    /* if y, convert output to Hertz */
	
    sf_hilbert_init(n1, nh, c);
    sf_deriv_init(n1, nh, c);
	
    mean=0.;
    for (i=i2=0; i2 < n2; i2++) {
	if (verb) sf_warning("slice %d of %d;",i2+1,n2);

	sf_floatread(trace,n1,in);
	sf_hilbert(trace,hilb);

	sf_deriv(trace,dtrace);
	sf_deriv(hilb,dhilb);

	sf_floatread(traced,n1,df);
	sf_hilbert(traced,hilbd);

	sf_deriv(traced,dtraced);
	sf_deriv(hilbd,dhilbd);

	for (i1=0; i1 < nh; i1++, i++) {
	    num[i] = 0.;
	    den[i] = 0.;
	    dnum[i] = 0.;
	    dden[i] = 0.;
	}	
	    
	for (i1=nh; i1 < n1-nh; i1++, i++) {
	    a = trace[i1];
	    b = hilb[i1];

	    da = traced[i1];
	    db = hilbd[i1];

	    num[i] = a*dhilb[i1]-b*dtrace[i1];
	    den[i] = a*a+b*b;

	    dnum[i] = da*dhilb[i1] + a*dhilbd[i1] - b*dtraced[i1] - db*dtrace[i1];
	    dden[i] = 2*a*da + 2*b*db;
	    
	    mean += den[i]*den[i];
	}
	    
	for (i1=n1-nh; i1 < n1; i1++, i++) {
	    num[i] = 0.;
	    den[i] = 0.;
	    dnum[i] = 0.;
	    dden[i] = 0.;
	}
    } /* i2 */
    
    if (verb) sf_warning(".");
    
    mean = sqrtf(n12/mean);
    
    for (i=0; i < n12; i++) {
	num[i] *= mean;
	den[i] *= mean;
    }

    sf_divn_init(dim, n12, n, rect, niter, true);
    sf_divn (num, den, phase);
    sf_divn_close();

    for (i=i2=0; i2 < n2; i2++) {
	if (verb) sf_warning("slice %d of %d;",i2+1,n2);

	for (i1=nh; i1 < n1-nh; i1++, i++) {
	    num[i] = 0.;
	}
	    
	for (i1=nh; i1 < n1-nh; i1++, i++) {
	    num[i] = (dnum[i] + dden[i]*phase[i])*mean;
	}

	for (i1=n1-nh; i1 < n1; i1++, i++) {
	    num[i] = 0.;
	}	
    } /* i2 */
    
    if (verb) sf_warning(".");

    sf_divn_init(dim, n12, n, rect, niter, true);
    sf_divn (num, den, phase);
    sf_divn_close();    
	
    if (hertz) {
	/* convert to Hertz */    
	d1 = 1./(2.*SF_PI*d1);
	for (i=0; i < n12; i++) {
	    phase[i] *= d1;
	}
    }
	
    sf_floatwrite(phase,n12,out);
	
    exit(0);
}


