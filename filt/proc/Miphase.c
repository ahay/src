/* Smooth estimate of instanteneous phase.
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

#include <rsf.h>

#include "divn.h"

int main (int argc, char* argv[])
{
    int n1, n2, i1, i2, i, n12, niter, dim, n[SF_MAX_DIM], rect[SF_MAX_DIM];
    float *trace, *num, *den, *phase, x1, x2, y1, y2;
    float complex *ctrace;
    char key[6];
    kiss_fft_cfg forw, invs;
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
    ctrace = sf_complexalloc (n1);

    forw = kiss_fft_alloc(n1,0,NULL,NULL);
    invs = kiss_fft_alloc(n1,1,NULL,NULL);
    if (NULL == forw || NULL == invs) sf_error("KISS FFT allocation error");

    num = sf_floatalloc(n12);
    den = sf_floatalloc(n12);
    phase = sf_floatalloc(n12);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(trace,n1,in);
	for (i1=0; i1 < n1; i1++) {
	    ctrace[i1]=trace[i1];
	}

	/* make an analytical phase */
	kiss_fft(forw,(const kiss_fft_cpx *) ctrace, (kiss_fft_cpx *) ctrace);
	ctrace[0] *= 0.5;
	ctrace[n1/2] *= 0.5;
	for (i1=n1/2+1; i1 < n1; i1++) {
	    ctrace[i1]=0.;
	}
	kiss_fft(invs,(const kiss_fft_cpx *) ctrace, (kiss_fft_cpx *) ctrace);

	x2 = y2 = 0.;
	for (i1=0; i1 < n1; i1++) {
	    x1 = x2;
	    y1 = y2;
	    x2 = crealf(ctrace[i1]); 
	    y2 = cimagf(ctrace[i1]); 

	    i = i1 + i2*n1;
	    num[i] = x1*y2-x2*y1;
	    den[i] = x1*x2+y1*y2;
	}
    }

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    divn_init(dim, n12, n, rect, niter);
    divn (num, den, phase);
    sf_floatwrite(phase,n12,out);

    exit(0);
}

/* 	$Id: Menvelope.c 696 2004-07-06 23:17:31Z fomels $	 */
