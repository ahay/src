/* Preconditioning for traveltime picking. 

Takes: rect1=1 rect2=1 ... 

rectN defines the size of the smoothing stencil in N-th dimension.
*/
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
    int shrt, lng;
    float *trace, *hilb, *num, *den, *phase, mean, c;
    char key[6];
    sf_triangle ts, tl;
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

    num = sf_floatalloc(n12);
    den = sf_floatalloc(n12);
    phase = sf_floatalloc(n12);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getint("order",&nh)) nh=10;
    /* Hilbert transformer order */
    if (!sf_getfloat("ref",&c)) c=1.;
    /* Hilbert transformer reference (0.5 < ref <= 1) */

    sf_hilbert_init(n1, nh, c);

    if (!sf_getint("short",&shrt)) shrt=1;
    /* short smoothing radius */
    if (!sf_getint("long",&lng)) lng=10;
    /* long smoothing radius */

    ts = sf_triangle_init(shrt,n1,false);
    tl = sf_triangle_init(lng,n1,false);

    mean=0.;
    for (i2=0; i2 < n2; i2++) {
	trace = num+n1*i2;
	hilb = den+n1*i2;

	sf_floatread(trace,n1,in);
	sf_hilbert(trace,hilb);

	for (i1=0; i1 < n1; i1++) {
	    trace[i1] = hypotf(trace[i1],hilb[i1]);
	    hilb[i1] = trace[i1];
	}

	sf_smooth2 (ts, 0, 1, false, trace);
	sf_smooth2 (tl, 0, 1, false, hilb);
	
	for (i1=0; i1 < nh; i1++) {
	    trace[i1] = 0.;
	    hilb[i1] = 0.;
	}
	for (i1=nh; i1 < n1-nh; i1++) {
	    mean += hilb[i1]*hilb[i1];
	}
	for (i1=n1-nh; i1 < n1; i1++) {
	    trace[i1] = 0.;
	    hilb[i1] = 0.;
	}
    }
    mean = sqrtf(n12/mean);

    for (i=0; i < n12; i++) {
	num[i] *= mean;
	den[i] *= mean;
    }

    sf_divn_init(dim, n12, n, rect, niter, true);
    sf_divn (num, den, phase);

    sf_floatwrite(phase,n12,out);

    exit(0);
}

/* 	$Id: Menvelope.c 696 2004-07-06 23:17:31Z fomels $	 */
