/* 3-D FFT encapsulated for wavefield extrapolation (OpenMP) */
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

#include "wexfft.h"

#include "wexutl.h"
/*^*/

/*------------------------------------------------------------*/
wexfft2d wexfft_init(wexcub3d cub,
		      int n1_, 
		      int n2_)
/*< initialize OMP fft >*/
{
    /*------------------------------------------------------------*/
    wexfft2d fft;
    fft = (wexfft2d) sf_alloc(1,sizeof(*fft));

    fft->n1 = n1_;
    fft->n2 = n2_;

    fft->ctmp1 = (kiss_fft_cpx*) sf_complexalloc(fft->n1);
    fft->ctmp2 = (kiss_fft_cpx*) sf_complexalloc(fft->n2);

    fft->forw1 = kiss_fft_alloc(fft->n1,0,NULL,NULL);
    fft->invs1 = kiss_fft_alloc(fft->n1,1,NULL,NULL);

    fft->forw2 = kiss_fft_alloc(fft->n2,0,NULL,NULL);
    fft->invs2 = kiss_fft_alloc(fft->n2,1,NULL,NULL);

    if (NULL == fft->forw2 || NULL == fft->invs2 || 
	NULL == fft->forw1 || NULL == fft->invs1) 
	sf_error("%s: KISS FFT allocation error",__FILE__);

    fft->fftscale = 1./sqrtf(fft->n1*fft->n2);

    return fft;
}

/*------------------------------------------------------------*/
void wexfft_close(wexfft2d fft)
/*< close OMP fft >*/
{
    free (fft->ctmp1);
    free (fft->ctmp2);

    free (fft->forw1);
    free (fft->invs1);

    free (fft->forw2);
    free (fft->invs2);
}

/*------------------------------------------------------------*/
void wexfft(bool inv          /* inverse/forward flag */, 
	     kiss_fft_cpx **pp /* [1...n2][1...n1] */,
	     wexfft2d fft) 
/*< apply 2-D FFT >*/
{
    int i1,i2;
    
    if (inv) {

	/* IFT 1 */
	for (i2=0; i2 < fft->n2; i2++) {
	    kiss_fft_stride(fft->invs1,
			    pp[i2],
			    fft->ctmp1,
			    1);

	    for (i1=0; i1<fft->n1; i1++) {
		pp[i2][i1] = fft->ctmp1[i1];
	    }
	}

	/* IFT 2 */
	for (i1=0; i1 < fft->n1; i1++) {
	    kiss_fft_stride(fft->invs2,
			    pp[0]+i1,
			    fft->ctmp2,
			    fft->n1);

	    for (i2=0; i2<fft->n2; i2++) {
		pp[i2][i1] = fft->ctmp2[i2];
	    }
	}

	/* scaling */
	for     (i2=0; i2<fft->n2; i2++) {
	    for (i1=0; i1<fft->n1; i1++) {
		pp[i2][i1] = sf_crmul(pp[i2][i1],fft->fftscale);
	    }
	}

    } else {

	/* scaling */
	for     (i2=0; i2<fft->n2; i2++) {
	    for (i1=0; i1<fft->n1; i1++) {
		pp[i2][i1] = sf_crmul(pp[i2][i1],fft->fftscale);
	    }
	}
	
	/* FFT 2 */
	for (i1=0; i1 < fft->n1; i1++) {
	    kiss_fft_stride(fft->forw2,
			    pp[0]+i1,
			    fft->ctmp2,
			    fft->n1);

	    for (i2=0; i2<fft->n2; i2++) {
		pp[i2][i1] = fft->ctmp2[i2];
	    }
	}

	/* FFT 1 */
	for (i2=0; i2 < fft->n2; i2++) {
	    kiss_fft_stride(fft->forw1,
			    pp[i2],
			    fft->ctmp1,
			    1);
	    
	    for (i1=0; i1<fft->n1; i1++) {
		pp[i2][i1] = fft->ctmp1[i1];
	    }
	}

    }

}

