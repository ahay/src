/* Complex 1-D FFT interface, complex to complex fft and ifft */
/*
  Copyright (C) 2010 University of Texas at Austin
 
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

static int n1;
static float wt;

static sf_complex *cc=NULL;
static kiss_fft_cfg cfg=NULL, icfg=NULL;

int cfft1_init(int nx             /* input data size */, 
	      int *nx2           /* padded data size */)
/*< initialize >*/
{
    int nk1;

    nk1 = n1 = kiss_fft_next_fast_size(nx);
    cc = sf_complexalloc(n1);
    *nx2 = n1;
    wt =  1.0/n1;

    return (nk1);
}

void cfft1(sf_complex *inp /* [n1] */, 
	  sf_complex *out /* [nk1] */)
/*< 1-D FFT >*/
{
    int i;

    if (NULL==cfg) {
      cfg  = kiss_fft_alloc(n1,0,NULL,NULL);
      if (NULL == cfg) sf_error("FFT allocation failure.");
    }

    for (i=0; i<n1; i++) {
#ifdef SF_HAS_COMPLEX_H
      cc[i] = i%2? inp[i]:(-1*inp[i]);
#else
      cc[i] = i%2? inp[i]:sf_cneg(inp[i]);
#endif
    }

    kiss_fft (cfg,(kiss_fft_cpx *) cc,(kiss_fft_cpx *) out);
    /*kiss_fft_stride(cfg,(kiss_fft_cpx *) cc,(kiss_fft_cpx *) out,1);*/

}

void icfft1(sf_complex *out /* [n1] */, 
  	    sf_complex *inp /* [nk1] */)
/*< 1-D inverse FFT >*/
{
    int i;

    if (NULL==icfg) {
      icfg = kiss_fft_alloc(n1,1,NULL,NULL);
      if (NULL == icfg) sf_error("FFT allocation failure.");
    }

    kiss_fft (icfg,(kiss_fft_cpx *) inp,(kiss_fft_cpx *) cc);
    /*kiss_fft_stride(icfg,(kiss_fft_cpx *) inp,(kiss_fft_cpx *) cc,1);*/

    for (i=0; i<n1; i++) {
#ifdef SF_HAS_COMPLEX_H
	out[i] = (i%2? wt:-wt) * cc[i];
#else
	out[i] = sf_crmul(cc[i],(i%2? wt:-wt));
#endif
    }
}
