/* Smooth estimate of instanteneous frequency. */
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

#include "divn.h"
#include "hilbert.h"

int main (int argc, char* argv[])
{
    int nh, n1,n2, i1,i2, i, n12, niter, dim, n[SF_MAX_DIM], rect[SF_MAX_DIM];
    float *trace, *hilb, *dtrace, *dhilb, *num, *den, *phase, a, b, c, mean;
    char key[6];
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    
    dim = sf_filedims (in,n);
    n1 = n[0];
    n12 = 1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	n12 *= n[i];
    }
    n2 = n12/n1;

    trace = sf_floatalloc(n1);
    hilb = sf_floatalloc(n1);
    dtrace = sf_floatalloc(n1);
    dhilb = sf_floatalloc(n1);

    num = sf_floatalloc(n12);
    den = sf_floatalloc(n12);
    phase = sf_floatalloc(n12);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getint("order",&nh)) nh=10;
    /* Hilbert transformer order */
    if (!sf_getfloat("ref",&c)) c=1.;
    /* Hilbert transformer reference (0.5 < ref <= 1) */

    hilbert_init(n1, nh, c);

    mean=0.;
    for (i=i2=0; i2 < n2; i2++) {
	sf_floatread(trace,n1,in);
	hilbert(trace,hilb);
	deriv(trace,dtrace);
	deriv(hilb,dhilb);

	for (i1=0; i1 < nh; i1++, i++) {
	    num[i] = 0.;
	    den[i] = 0.;
	}
	for (i1=nh; i1 < n1-nh; i1++, i++) {
	    a = trace[i1];
	    b = hilb[i1];
	    num[i] = a*dhilb[i1]-b*dtrace[i1];
	    den[i] = a*a+b*b;
	    mean += num[i];
	}
	for (i1=n1-nh; i1 < n1; i1++, i++) {
	    num[i] = 0.;
	    den[i] = 0.;
	}
    }
    mean = sqrtf(n12/mean);

    for (i=0; i < n12; i++) {
	num[i] *= mean;
	den[i] *= mean;
    }

    divn_init(dim, n12, n, rect, niter);
    divn (num, den, phase);
    sf_floatwrite(phase,n12,out);

    exit(0);
}

/* 	$Id: Menvelope.c 696 2004-07-06 23:17:31Z fomels $	 */
