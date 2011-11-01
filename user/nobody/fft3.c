/* 3-D FFT encapsulated */
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

#include "fft3.h"

static kiss_fft_cfg forw1, invs1; /* FFT on axis 1 */
static kiss_fft_cfg forw2, invs2; /* FFT on axis 2 */
static kiss_fft_cfg forw3, invs3; /* FFT on axis 3 */
static int          n1,n2,n3;
static float        fftscale;
static kiss_fft_cpx *trace2, *trace3;
static sf_complex   *shf1,*shf2,*shf3;

/*------------------------------------------------------------*/
void fft3_init(int n1_, int n2_, int n3_)
/*< initialize >*/
{
    n1 = n1_; 
    n2 = n2_;
    n3 = n3_;

    forw1 = kiss_fft_alloc(n1,0,NULL,NULL);
    invs1 = kiss_fft_alloc(n1,1,NULL,NULL);

    forw2 = kiss_fft_alloc(n2,0,NULL,NULL);
    invs2 = kiss_fft_alloc(n2,1,NULL,NULL);

    forw3 = kiss_fft_alloc(n3,0,NULL,NULL);
    invs3 = kiss_fft_alloc(n3,1,NULL,NULL);
    
    if (NULL == forw1 || NULL == invs1 || 
	NULL == forw2 || NULL == invs2 ||
	NULL == forw3 || NULL == invs3) 
	sf_error("%s: KISS FFT allocation error",__FILE__);

    trace2 = (kiss_fft_cpx*) sf_complexalloc(n2);
    trace3 = (kiss_fft_cpx*) sf_complexalloc(n3);

    fftscale = 1./sqrtf(n1*n2*n3);
}

/*------------------------------------------------------------*/
void fft3_close(void)
/*< Free allocated storage >*/
{
    free (trace2);
    free (trace3);

    free (forw1);
    free (invs1);
    free (forw2);
    free (invs2);
    free (forw3);
    free (invs3);
}

/*------------------------------------------------------------*/
void fft3(bool inv           /* inverse/forward flag */, 
	  kiss_fft_cpx ***pp /* [n1][n2][n3] */) 
/*< Apply 3-D FFT >*/
{
    int i1, i2, i3;
    
    if (inv) {

	/* IFT 1 */
	for    (i3=0; i3 < n3; i3++) {
	    for(i2=0; i2 < n2; i2++) {
		kiss_fft(invs1,pp[i3][i2],pp[i3][i2]);
	    }
	}

	/* IFT 2 */
	for    (i3=0; i3 < n3; i3++) {
	    for(i1=0; i1 < n1; i1++) {
		kiss_fft_stride(invs2,pp[i3][0]+i1,trace2,n1);
		for(i2=0; i2 < n2; i2++) {
		    pp[i3][i2][i1] = trace2[i2];
		}
	    }
	}
	
	/* IFT 3 */
	for    (i2=0; i2 < n2; i2++) {
	    for(i1=0; i1 < n1; i1++) {
		kiss_fft_stride(invs3,pp[0][0]+i1+i2*n1,trace3,n1*n2);
		for(i3=0; i3 < n3; i3++) {
		    pp[i3][i2][i1] = trace3[i3];
		}
	    }
	}

	/* scaling */
	for        (i3=0; i3 < n3; i3++) {
	    for    (i2=0; i2 < n2; i2++) {
		for(i1=0; i1 < n1; i1++) {
		    pp[i3][i2][i1] = sf_crmul(pp[i3][i2][i1],fftscale);
		}
	    }
	}
	
    } else {
	/* scaling */
	for        (i3=0; i3 < n3; i3++) {
	    for    (i2=0; i2 < n2; i2++) {
		for(i1=0; i1 < n1; i1++) {
		    pp[i3][i2][i1] = sf_crmul(pp[i3][i2][i1],fftscale);
		}
	    }
	}

	/* FFT 3 */
	for    (i2=0; i2 < n2; i2++) {
	    for(i1=0; i1 < n1; i1++) {
		kiss_fft_stride(forw3,pp[0][0]+i1+i2*n1,trace3,n1*n2);
		for(i3=0; i3 < n3; i3++) {
		    pp[i3][i2][i1] = trace3[i3];
		}
	    }
	}

	/* FFT 2 */
	for    (i3=0; i3 < n3; i3++) {
	    for(i1=0; i1 < n1; i1++) {
		kiss_fft_stride(forw2,pp[i3][0]+i1,trace2,n1);
		for(i2=0; i2 < n2; i2++) {
		    pp[i3][i2][i1] = trace2[i2];
		}
	    }
	}

	/* FFT 1 */
	for    (i3=0; i3 < n3; i3++) {
	    for(i2=0; i2 < n2; i2++) {
		kiss_fft(forw1,pp[i3][i2],pp[i3][i2]);
	    }
	}

    }
}

/*------------------------------------------------------------*/
void sft3_init(float o1, float d1, 
	       float o2, float d2,
	       float o3, float d3)
/*< origin shift (assumes no centering) >*/
{
    int   i1,i2,i3;
    int   k1,k2,k3;
    float w1,w2,w3;
    float shift;

    w1=2.0*SF_PI/(n1*d1) * o1;
    w2=2.0*SF_PI/(n2*d2) * o2;
    w3=2.0*SF_PI/(n3*d3) * o3;

    k1=n1/2;
    k2=n2/2;
    k3=n3/2;

    shf1 = sf_complexalloc(n1);
    for(i1=0; i1<n1; i1++) { shf1[i1]=1.0; }

    for(i1=0; i1<k1; i1++) {
	shift = w1 * i1;
	shf1[i1]    = sf_cmplx(cosf(shift),sinf(shift));

	shift = w1 * (-k1-1+i1);
	shf1[k1+i1] = sf_cmplx(cosf(shift),sinf(shift));
    }

    shf2 = sf_complexalloc(n2);
    for(i2=0; i2<n2; i2++) { shf2[i2]=1.0; }

    for(i2=0; i2<k2; i2++) {
	shift = w2 * i2;
	shf2[i2]    = sf_cmplx(cosf(shift),sinf(shift));

	shift = w2 * (-k2-1+i2);
	shf2[k2+i2] = sf_cmplx(cosf(shift),sinf(shift));
    }

    shf3 = sf_complexalloc(n3);
    for(i3=0; i3<n3; i3++) { shf3[i3]=1.0; }

    for(i3=0; i3<k3; i3++) {
	shift = w3 * i3;
	shf3[i3]    = sf_cmplx(cosf(shift),sinf(shift));

	shift = w3 * (-k3-1+i3);
	shf2[k3+i3] = sf_cmplx(cosf(shift),sinf(shift));
    }
}

/*------------------------------------------------------------*/
void sft3_close()
/*< close shift >*/
{
    free(shf1);
    free(shf2);
    free(shf3);
}

/*------------------------------------------------------------*/
void sft3(sf_complex ***pp)
/*< apply shift >*/
{
    int i1,i2,i3;

    for        (i3=0; i3<n3; i3++) {
	for    (i2=0; i2<n2; i2++) {
	    for(i1=0; i1<n1; i1++) {
#ifdef SF_HAS_COMPLEX_H
		pp[i3][i2][i1] *= shf1[i1]*shf2[i2]*shf3[i3];
#else
		pp[i3][i2][i1] = sf_cmul(pp[i3][i2][i1],
					 sf_cmul(sf_cmul(shf1[i1],shf2[i2]),shf3[i3])
		    );
#endif
	    }
	}    
    }

}
