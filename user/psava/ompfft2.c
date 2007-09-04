/* 2-D FFT encapsulated (OpenMP version) */
/*
  Copyright (C) 2007 Colorado School of Mines
  
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

#include "ompfft2.h"

static kiss_fft_cfg *forw1;
static kiss_fft_cfg *invs1; /* FFT on axis 1 */
static kiss_fft_cfg *forw2;
static kiss_fft_cfg *invs2; /* FFT on axis 2 */

static kiss_fft_cpx **ctrace;

static int          n1,n2;
static float        fftscale;
static sf_complex *shf1,*shf2;

/*------------------------------------------------------------*/
void ompfft2_init(int n1_, 
		  int n2_,
		  int ompnth)
/*< initialize >*/
{
    int ompith;
    n1 = n1_;
    n2 = n2_;

    ctrace = (kiss_fft_cpx**) sf_complexalloc2(n2,ompnth);

    forw1 = (kiss_fft_cfg*) sf_alloc(ompnth,sizeof(kiss_fft_cfg));
    invs1 = (kiss_fft_cfg*) sf_alloc(ompnth,sizeof(kiss_fft_cfg));
    forw2 = (kiss_fft_cfg*) sf_alloc(ompnth,sizeof(kiss_fft_cfg));
    invs2 = (kiss_fft_cfg*) sf_alloc(ompnth,sizeof(kiss_fft_cfg));

    for(ompith=0; ompith<ompnth; ompith++) {
	forw1[ompith] = kiss_fft_alloc(n1,0,NULL,NULL);
	invs1[ompith] = kiss_fft_alloc(n1,1,NULL,NULL);
	forw2[ompith] = kiss_fft_alloc(n2,0,NULL,NULL);
	invs2[ompith] = kiss_fft_alloc(n2,1,NULL,NULL);

	if (NULL == forw2[ompith] || NULL == invs2[ompith] || 
	    NULL == forw1[ompith] || NULL == invs1[ompith]) 
	    sf_error("%s: KISS FFT allocation error",__FILE__);
    }

    fftscale = 1./sqrtf(n1*n2);
}

/*------------------------------------------------------------*/
void ompfft2_close(void)
/*< Free allocated storage >*/
{
    free(*ctrace); free (ctrace);

    free (*forw2); free (forw2);
    free (*invs2); free (invs2);
    free (*forw1); free (forw1);
    free (*invs1); free (invs1);
}

/*------------------------------------------------------------*/
void ompfft2(bool inv          /* inverse/forward flag */, 
	     kiss_fft_cpx **pp /* [1...n2][1...n1] */,
	     int ompith) 
/*< Apply 2-D FFT >*/
{
    int i1,i2;
    
    if (inv) {
	for (i2=0; i2 < n2; i2++) {
#pragma omp critical
	    kiss_fft(invs1[ompith],pp[i2],pp[i2]);
	}
	for (i1=0; i1 < n1; i1++) {
#pragma omp critical
	    kiss_fft_stride(invs2[ompith],pp[0]+i1,ctrace[ompith],n1);
	    for (i2=0; i2<n2; i2++) {
		pp[i2][i1] = ctrace[ompith][i2];
	    }
	}
	
	for (i2=0; i2<n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		pp[i2][i1] = sf_crmul(pp[i2][i1],fftscale);
	    }
	}
    } else {
	for (i2=0; i2<n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		pp[i2][i1] = sf_crmul(pp[i2][i1],fftscale);
	    }
	}
	
	for (i1=0; i1 < n1; i1++) {
#pragma omp critical
	    kiss_fft_stride(forw2[ompith],pp[0]+i1,ctrace[ompith],n1);
	    for (i2=0; i2<n2; i2++) {
		pp[i2][i1] = ctrace[ompith][i2];
	    }
	}
	for (i2=0; i2 < n2; i2++) {
#pragma omp critical
	    kiss_fft(forw1[ompith],pp[i2],pp[i2]);
	}
    }
}

/*------------------------------------------------------------*/
void ompsft2_init(float o1, float d1, float o2, float d2)
/*< origin shift >*/
{
    int i1,i2;
    float shift;

    shf1 = sf_complexalloc(n1);
    for( i1=0; i1<n1; i1++) {
	shift = 2.0*SF_PI*i1/n1*o1/d1;
	shf1[i1] = sf_cmplx(cosf(shift),sinf(shift));
    }

    shf2 = sf_complexalloc(n2);
    for( i2=0; i2<n2; i2++) {
	shift = 2.0*SF_PI*i2/n2*o2/d2;
	shf2[i2] = sf_cmplx(cosf(shift),sinf(shift));
    }
}

/*------------------------------------------------------------*/
void ompsft2_close()
/*< close shift >*/
{
    free(shf1);
    free(shf2);
}

/*------------------------------------------------------------*/
void ompsft2(sf_complex **pp)
/*< apply shift >*/
{
    int i1,i2;

    for(i2=0; i2<n2; i2++) {
	for (i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
	    pp[i2][i1] *= shf1[i1]*shf2[i2];
#else
	    pp[i2][i1] = sf_cmul(pp[i2][i1],sf_cmul(shf1[i1],shf2[i2]));
#endif
	}
    }    
}

/*------------------------------------------------------------*/
void ompcnt2(sf_complex **pp)
/*< apply centering >*/
{
    int i1,i2;

    for(i2=1; i2<n2; i2+=2) {
	for(i1=1; i1<n1; i1+=2) {
#ifdef SF_HAS_COMPLEX_H
	    pp[i2][i1] = - pp[i2][i1];
#else
	    pp[i2][i1] = sf_cneg(pp[i2][i1]);
#endif
	}
    }
}
