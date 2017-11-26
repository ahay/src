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

#ifdef SF_HAS_FFTW
#include <fftw3.h>
#endif

#include "alloc.h"
#include "kiss_fftr.h"
#include "cosft.h"
#include "komplex.h"

static int nt, nw, n1;
static float *p /* , dt */;
static kiss_fft_cpx *pp;

#ifdef SF_HAS_FFTW
static fftwf_plan cfg, icfg;
#else
static kiss_fftr_cfg forw, invs;
#endif

void sf_cosft_init(int n1_in)
/*< initialize >*/ 
{
    n1 = n1_in;
    nt = 2*kiss_fft_next_fast_size(n1-1);
    nw = nt/2+1;
    p  = sf_floatalloc (nt);
    pp = (kiss_fft_cpx*) sf_complexalloc(nw);

#ifdef SF_HAS_FFTW
    cfg = fftwf_plan_dft_r2c_1d(nt, p, (fftwf_complex *) pp,
				FFTW_ESTIMATE);
    icfg = fftwf_plan_dft_c2r_1d(nt, (fftwf_complex *) pp, p,
				 FFTW_ESTIMATE);
#else
    forw = kiss_fftr_alloc(nt,0,NULL,NULL);
    invs = kiss_fftr_alloc(nt,1,NULL,NULL);
#endif
}

void sf_cosft_close(void) 
/*< free allocated storage >*/
{
    free (p);
    free (pp);
#ifdef SF_HAS_FFTW
    fftwf_destroy_plan(cfg);
    fftwf_destroy_plan(icfg);
#else
    free (forw);
    free (invs);
#endif
}

void sf_cosft_frw (float *q /* data */, 
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

#ifdef SF_HAS_FFTW   
    fftwf_execute(cfg);
#else
    kiss_fftr(forw, p, pp);
#endif
    
    for (i=0; i < n1; i++) {
	q[o1+i*d1] = sf_crealf(pp[i]);
    }
}

void sf_cosft_inv (float *q /* data */, 
		   int o1   /* first sample */, 
		   int d1   /* step */) 
/*< inverse transform >*/
{
    int i;
	
	
    for (i=0; i < n1; i++) {
	pp[i].r = q[o1+i*d1];
	pp[i].i = 0.;
    }
    /* pad */
    for (i=n1; i < nw; i++) { 
	pp[i].r = 0.; 
	pp[i].i = 0.;
    }

#ifdef SF_HAS_FFTW   
    fftwf_execute(icfg);
#else
    kiss_fftri(invs,pp,p);
#endif
    
    for (i=0; i < n1; i++) {
	q[o1+i*d1] = p[i]/nt;
    }
}

/* 	$Id$	 */
