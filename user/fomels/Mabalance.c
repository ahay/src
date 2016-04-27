/* Amplitude balancing. */
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
    bool reverse;
    int nh, n1,n2, i1,i2, i, n12, niter, dim, n[SF_MAX_DIM], rect[SF_MAX_DIM];
    float *trace, *hilb, *trace2, *hilb2, *num, *den, *rat, *org, c, mean;
    char key[6];
    sf_file in, out, ref, weight;

    sf_init (argc,argv);
    in = sf_input("in");
    ref = sf_input("other");
    out = sf_output("out");

    if (NULL != sf_getstring("weight")) { /* optional weight output */
	weight = sf_output("weight");
    } else {
	weight = NULL;
    }

    if (!sf_getbool("reverse",&reverse)) reverse=true;
    /* reverse weight */

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
    trace2 = sf_floatalloc(n1);
    hilb2 = sf_floatalloc(n1);

    org = sf_floatalloc(n12);
    num = sf_floatalloc(n12);
    den = sf_floatalloc(n12);
    rat = sf_floatalloc(n12);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getint("order",&nh)) nh=100;
    /* Hilbert transformer order */
    if (!sf_getfloat("ref",&c)) c=1.;
    /* Hilbert transformer reference (0.5 < ref <= 1) */

    sf_hilbert_init(n1, nh, c);

    mean=0.;
    for (i=i2=0; i2 < n2; i2++) {
	sf_floatread(trace,n1,in);
	sf_hilbert(trace,hilb);

	sf_floatread(trace2,n1,ref);
	sf_hilbert(trace2,hilb2);

	for (i1=0; i1 < nh; i1++, i++) {
	    num[i] = 0.;
	    den[i] = 0.;
	    org[i] = trace[i1];
	}
	for (i1=nh; i1 < n1-nh; i1++, i++) {
	    if (reverse) {
		num[i] = hypotf(trace[i1],hilb[i1]);
		den[i] = hypotf(trace2[i1],hilb2[i1]);
	    } else {
		num[i] = hypotf(trace2[i1],hilb2[i1]);
		den[i] = hypotf(trace[i1],hilb[i1]);
	    }
	    org[i] = trace[i1];
	    mean += den[i];
	}
	for (i1=n1-nh; i1 < n1; i1++, i++) {
	    num[i] = 0.;
	    den[i] = 0.;
	    org[i] = trace[i1];
	}
    }
    mean = sqrtf(n12/mean);

    for (i=0; i < n12; i++) {
	num[i] *= mean;
	den[i] *= mean;
    }

    sf_divn_init(dim, n12, n, rect, niter, true);
    sf_divn (num, den, rat);

    for (i=0; i < n12; i++) {
	if (reverse) {
	    if (rat[i] != 0.0f) {
		rat[i] = 1.0f/rat[i];
		org[i] *= rat[i];
	    }
	} else {
	    org[i] *= rat[i];
	}
    }
    
    sf_floatwrite(org,n12,out);
    if (NULL != weight) sf_floatwrite(rat,n12,weight);

    exit(0);
}

