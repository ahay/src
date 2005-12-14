/* Cosine Fourier transform */
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

#include "cosft.h"

static int nt, nw, n1;
static float *p /* , dt */;
static float complex *pp;
static kiss_fftr_cfg forw, invs;

void cosft_init(int n1_in)
/*< initialize >*/ 
{
    n1 = n1_in;
    nt = sf_fftr_size(2*(n1-1));
    nw = nt/2+1;
    p  = sf_floatalloc (nt);
    pp = sf_complexalloc(nw);
    forw = kiss_fftr_alloc(nt,0,NULL,NULL);
    invs = kiss_fftr_alloc(nt,1,NULL,NULL);
}

void cosft_close(void) 
/*< free allocated storage >*/
{
    free (p);
    free (pp);
    free (forw);
    free (invs);
}

void cosft_frw (float *q /* data */, 
		int o1   /* first sample */, 
		int d1   /* step */) 
/*< forward transform >*/
{
    int i;

    for (i=0; i < n1; i++) {
	p[i] = q[o1+i*d1];
    }
    for (i=n1; i < nw; i++) { 
	p[i] = 0.; /* pad */
    }
    for (i=nw; i < nt; i++) {
	p[i] = p[nt-i];
    }
    
    kiss_fftr(forw, p, (kiss_fft_cpx *) pp);
    
    for (i=0; i < n1; i++) {
	q[o1+i*d1] = crealf(pp[i]);
    }
}

void cosft_inv (float *q /* data */, 
		int o1   /* first sample */, 
		int d1   /* step */) 
/*< inverse transform >*/
{
    int i;

    for (i=0; i < n1; i++) {
	pp[i] = q[o1+i*d1];
    }
    for (i=n1; i < nw; i++) { 
	pp[i] = 0.; /* pad */
    }

/*

    if (0. != dt) {
	for (i=0; i < n; i++) {
	    pp[i] *= cexpf(I*i*dt);
	}
    }
  
*/
  
    kiss_fftri(invs,(const kiss_fft_cpx *) pp, p);
    
    for (i=0; i < n1; i++) {
	q[o1+i*d1] = p[i]/nt;
    }
}

/* 	$Id$	 */
