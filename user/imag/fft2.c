/* 2-D FFT encapsulated */
/*
  Copyright (C) 2008 Colorado School of Mines
  
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
static int          n1,n2;
static float        fftscale;
static kiss_fft_cpx *trace2;
static sf_complex   *shf1,*shf2;

/*------------------------------------------------------------*/
void fft2_init(int n1_, int n2_)
/*< initialize >*/
{
    n1 = n1_;
    n2 = n2_;

    forw1 = kiss_fft_alloc(n1,0,NULL,NULL);
    invs1 = kiss_fft_alloc(n1,1,NULL,NULL);

    forw2 = kiss_fft_alloc(n2,0,NULL,NULL);
    invs2 = kiss_fft_alloc(n2,1,NULL,NULL);

    trace2 = (kiss_fft_cpx*) sf_complexalloc(n2);

    if (NULL == forw2 || NULL == invs2 || 
	NULL == forw1 || NULL == invs1) 
	sf_error("%s: KISS FFT allocation error",__FILE__);

    fftscale = 1./sqrtf(n1*n2);
}

/*------------------------------------------------------------*/
void fft2_close(void)
/*< Free allocated storage >*/
{
    free (trace2);

    free (forw1);
    free (invs1);
    free (forw2);
    free (invs2);
}

/*------------------------------------------------------------*/
void fft2(bool inv          /* inverse/forward flag */, 
	  kiss_fft_cpx **pp /* [1...n2][1...n1] */) 
/*< Apply 2-D FFT >*/
{
    int i1,i2;
    
    if (inv) {

	/* IFT 1 */
	for(i2=0; i2 < n2; i2++) {
	    kiss_fft(invs1,pp[i2],pp[i2]);
	}

	/* IFT 2 */
	for(i1=0; i1 < n1; i1++) {
	    kiss_fft_stride(invs2,pp[0]+i1,trace2,n1);
	    for(i2=0; i2<n2; i2++) {
		pp[i2][i1] = trace2[i2];
	    }
	}

	/* scaling */
	for    (i2=0; i2<n2; i2++) {
	    for(i1=0; i1 < n1; i1++) {
		pp[i2][i1] = sf_crmul(pp[i2][i1],fftscale);
	    }
	}
    } else {

	/* scaling */
	for    (i2=0; i2<n2; i2++) {
	    for(i1=0; i1 < n1; i1++) {
		pp[i2][i1] = sf_crmul(pp[i2][i1],fftscale);
	    }
	}

	/* FFT 2 */
	for(i1=0; i1 < n1; i1++) {
	    kiss_fft_stride(forw2,pp[0]+i1,trace2,n1);
	    for(i2=0; i2<n2; i2++) {
		pp[i2][i1] = trace2[i2];
	    }
	}

	/* FFT 1 */
	for(i2=0; i2 < n2; i2++) {
	    kiss_fft(forw1,pp[i2],pp[i2]);
	}
    }
}

/*------------------------------------------------------------*/
void sft2_init(float o1, float d1, 
	       float o2, float d2)
/*< origin shift (assumes no centering) >*/
{
    int   i1,i2;
    int   k1,k2;
    float w1,w2;
    float shift;

    w1=2.0*SF_PI/(n1*d1) * o1;
    w2=2.0*SF_PI/(n2*d2) * o2;

    k1=n1/2;
    k2=n2/2;

    shf1 = sf_complexalloc(n1);
    for(i1=0; i1<n1; i1++) { shf1[i1]=sf_cmplx(1.0,0.0); }

    for(i1=0; i1<k1; i1++) {
	shift = w1 * i1;
	shf1[i1]    = sf_cmplx(cosf(shift),sinf(shift));

	shift = w1 * (-k1-1+i1);
	shf1[k1+i1] = sf_cmplx(cosf(shift),sinf(shift));
    }

    shf2 = sf_complexalloc(n2);
    for(i2=0; i2<n2; i2++) { shf2[i2]=sf_cmplx(1.0,0.0); }

    for(i2=0; i2<k2; i2++) {
	shift = w2 * i2;
	shf2[i2]    = sf_cmplx(cosf(shift),sinf(shift));

	shift = w2 * (-k2-1+i2);
	shf2[k2+i2] = sf_cmplx(cosf(shift),sinf(shift));
    }
}

/*------------------------------------------------------------*/
void sft2_close()
/*< close shift >*/
{
    free(shf1);
    free(shf2);
}

/*------------------------------------------------------------*/
void sft2(sf_complex **pp)
/*< apply shift >*/
{
    int i1,i2;

    for    (i2=0; i2<n2; i2++) {
	for(i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
	    pp[i2][i1] *= shf1[i1]*shf2[i2];
#else
	    pp[i2][i1] = sf_cmul(pp[i2][i1],
				 sf_cmul(shf1[i1],shf2[i2])
		);
#endif
	}
    }    
}

/*------------------------------------------------------------*/
void cnt2(sf_complex **pp)
/*< apply centering >*/
{
    int i1,i2;

    for    (i2=1; i2<n2; i2+=2) {
	for(i1=1; i1<n1; i1+=2) {
#ifdef SF_HAS_COMPLEX_H
	    pp[i2][i1] = - pp[i2][i1];
#else
	    pp[i2][i1] = sf_cneg(pp[i2][i1]);
#endif
	}
    }
}
