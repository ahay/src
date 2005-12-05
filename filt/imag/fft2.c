/* 2-D FFT encapsulated */
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
/*^*/

#include "fft2.h"

static kiss_fft_cfg forw1, invs1; /* FFT on axis 1 */
static kiss_fft_cfg forw2, invs2; /* FFT on axis 2 */
static int            n1,n2;
static float          fftscale;
static float complex *ctrace;
static float complex *shf1,*shf2;

void fft2_init(int n1_, int n2_)
/*< initialize >*/
{
    n1 = n1_;
    n2 = n2_;

    forw1 = kiss_fft_alloc(n1,0,NULL,NULL);
    invs1 = kiss_fft_alloc(n1,1,NULL,NULL);
    forw2 = kiss_fft_alloc(n2,0,NULL,NULL);
    invs2 = kiss_fft_alloc(n2,1,NULL,NULL);

    ctrace = sf_complexalloc(n2);

    if (NULL == forw2 || NULL == invs2 || 
	NULL == forw1 || NULL == invs1) 
	sf_error("%s: KISS FFT allocation error",__FILE__);

    fftscale = 1./sqrtf(n1*n2);
}

void fft2_close(void)
/*< Free allocated storage >*/
{
    free (ctrace);
    free (forw2);
    free (invs2);
    free (forw1);
    free (invs1);
}

void fft2(bool inv           /* inverse/forward flag */, 
	  complex float **pp /* [1...n2][1...n1] */) 
/*< Apply 2-D FFT >*/
{
    int i1,i2;
    
    if (inv) {
	for (i2=0; i2 < n2; i2++) {
	    kiss_fft(invs1,
		     (const kiss_fft_cpx *) pp[i2], 
		     (      kiss_fft_cpx *) pp[i2]);
	}
	for (i1=0; i1 < n1; i1++) {
	    kiss_fft_stride(invs2,
			    (const kiss_fft_cpx *) (pp[0]+i1), 
			    (      kiss_fft_cpx *) ctrace,n1);
	    for (i2=0; i2<n2; i2++) {
		pp[i2][i1] = ctrace[i2];
	    }
	}

	for (i2=0; i2<n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		pp[i2][i1] *= fftscale;
	    }
	}
    } else {
	for (i2=0; i2<n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		pp[i2][i1] *= fftscale;
	    }
	}

	for (i1=0; i1 < n1; i1++) {
	    kiss_fft_stride(forw2,
			    (const kiss_fft_cpx *) (pp[0]+i1), 
			    (      kiss_fft_cpx *) ctrace,n1);
	    for (i2=0; i2<n2; i2++) {
		pp[i2][i1] = ctrace[i2];
	    }
	}
	for (i2=0; i2 < n2; i2++) {
	    kiss_fft(forw1,
		     (const kiss_fft_cpx *) pp[i2], 
		     (      kiss_fft_cpx *) pp[i2]);
	}
    }
}

void sft2_init(float o1, float d1, float o2, float d2)
/*< origin shift >*/
{
    int i1,i2;

    shf1 = sf_complexalloc(n1);
    for( i1=0; i1<n1; i1++) {
	shf1[i1] = cexpf(+I*2.0*SF_PI*i1/n1*o1/d1);
    }

    shf2 = sf_complexalloc(n2);
    for( i2=0; i2<n2; i2++) {
	shf2[i2] = cexpf(+I*2.0*SF_PI*i2/n2*o2/d2);
    }
}

void sft2_close()
/*< close shift >*/
{
    free(shf1);
    free(shf2);
}

void sft2(complex float **pp)
/*< apply shift >*/
{
    int i1,i2;

    for(i2=0; i2<n2; i2++) {
	for (i1=0; i1<n1; i1++) {
	    pp[i2][i1] *= shf1[i1]*shf2[i2];
	}
    }    
}

void cnt2(complex float **pp)
/*< apply centering >*/
{
    int i1,i2;

    for(i2=1; i2<n2; i2+=2) {
	for(i1=1; i1<n1; i1+=2) {
	    pp[i2][i1] = - pp[i2][i1];
	}
    }
}
