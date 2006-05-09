/* 1-D FFT encapsulated */
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

#include "fft1.h"

static kiss_fft_cfg forw1, invs1; /* FFT on axis 1 */
static kiss_fft_cfg forw2, invs2; /* FFT on axis 2 */
static int            n1,n2;
static int            axis;
static float          fftscale;
static kiss_fft_cpx *ctrace;
static sf_complex *shf1,*shf2;

void fft1_init(int n1_, int n2_, int axis_)
/*< initialize >*/
{
    n1 = n1_;
    n2 = n2_;
    axis = axis_;

    if(axis==1) {
	forw1 = kiss_fft_alloc(n1,0,NULL,NULL);
	invs1 = kiss_fft_alloc(n1,1,NULL,NULL);

	if (NULL == forw1 || NULL == invs1) 
	    sf_error("%s: KISS FFT allocation error (axis 1)",__FILE__);

	fftscale = 1./n1;
    } else {
	ctrace = (kiss_fft_cpx*) sf_complexalloc(n2);

	forw2 = kiss_fft_alloc(n2,0,NULL,NULL);
	invs2 = kiss_fft_alloc(n2,1,NULL,NULL);

	if (NULL == forw2 || NULL == invs2) 
	    sf_error("%s: KISS FFT allocation error (axis 2)",__FILE__);

	fftscale = 1./n2;
    }
}

void fft1_close(void)
/*< Free allocated storage >*/
{
    if(axis==1) {
	free (forw1);
	free (invs1);
    } else {
	free (forw2);
	free (invs2);
	free (ctrace);
    }
}

void fft1a1(bool inv, 
	    kiss_fft_cpx **pp /* [1...n2][1...n1] */) 
/*< Apply 1-D FFT on axis 1 >*/
{
    int i1,i2;
    
    if (inv) {
	for (i2=0; i2 < n2; i2++) {
	    kiss_fft(invs1,pp[i2],pp[i2]);
	}
	for (i1=0; i1 < n1; i1++) {
	    for (i2=0; i2<n2; i2++) {
		pp[i2][i1] = sf_crmul(pp[i2][i1],fftscale);
	    }
	}
    } else {
	for (i2=0; i2 < n2; i2++) {
	    kiss_fft(forw1,pp[i2],pp[i2]);
	}
    }
}

void fft1a2(bool inv, 
	    kiss_fft_cpx **pp /* [1...n2][1...n1] */) 
/*< Apply 1-D FFT on axis 2 >*/
{
    int i1,i2;
    
    if (inv) {
	for (i1=0; i1 < n1; i1++) {
	    kiss_fft_stride(invs2,pp[0]+i1,ctrace,n1);
	    for (i2=0; i2<n2; i2++) {
		pp[i2][i1] = ctrace[i2];
	    }
	}
	for (i1=0; i1 < n1; i1++) {
	    for (i2=0; i2<n2; i2++) {
		pp[i2][i1] = sf_crmul(pp[i2][i1],fftscale);
	    }
	}
    } else {
	for (i1=0; i1 < n1; i1++) {
	    kiss_fft_stride(forw2,pp[0]+i1,ctrace,n1);
	    for (i2=0; i2<n2; i2++) {
		pp[i2][i1] = ctrace[i2];
	    }
	}
    }
}

/*------------------------------------------------------------*/

void sft1_init(float o, float d)
/*< origin shift >*/
{
    int i1,i2;
    float t;

    if(axis==1) {
	shf1 = sf_complexalloc(n1);
	for( i1=0; i1<n1; i1++) {
	    t = 2.0*SF_PI*i1/n1*o/d;
	    shf1[i1] = sf_cmplx(cosf(t),sinf(t));
	}
    } else {
	shf2 = sf_complexalloc(n2);
	for( i2=0; i2<n2; i2++) {
	    t = 2.0*SF_PI*i2/n2*o/d;
	    shf2[i2] = sf_cmplx(cosf(t),sinf(t));
	}
    }
}

void sft1_close()
/*< close shift >*/
{
    if(axis==1) {
	free(shf1);
    } else {
	free(shf2);
    }
}

void sft1a1(sf_complex **pp)
/*< apply shift >*/
{
    int i1,i2;

    for(i2=0; i2<n2; i2++) {
	for (i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
	    pp[i2][i1] *= shf1[i1];
#else
	    pp[i2][i1] = sf_cmul(pp[i2][i1],shf1[i1]);
#endif
	}
    }
}

void sft1a2(sf_complex **pp)
/*< apply shift >*/
{
    int i1,i2;

    for(i2=0; i2<n2; i2++) {
	for(i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
	    pp[i2][i1] *= shf2[i2];
#else
	    pp[i2][i1] = sf_cmul(pp[i2][i1],shf2[i2]);
#endif
	}
    }
}

/*------------------------------------------------------------*/

void cnt1a1(sf_complex **pp)
/*< apply centering >*/
{
    int i1,i2;

    for(i2=0; i2<n2; i2++) {
	for(i1=1; i1<n1; i1+=2) {
#ifdef SF_HAS_COMPLEX_H
	    pp[i2][i1] = - pp[i2][i1];
#else
	    pp[i2][i1] = sf_cneg(pp[i2][i1]);
#endif
	}
    }
}

void cnt1a2(sf_complex **pp)
/*< apply centering >*/
{
    int i1,i2;

    for(i2=1; i2<n2; i2+=2) {
	for(i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
	    pp[i2][i1] = - pp[i2][i1];
#else
	    pp[i2][i1] = sf_cneg(pp[i2][i1]);
#endif
	}
    }
}
