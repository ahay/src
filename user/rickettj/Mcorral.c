/* Cross-correlate every trace with every other in frequency domain. */
/*
  Copyright (C) 2000 The Board of Trustees of Stanford University
  
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

#ifdef SF_HAS_FFTW
#include <fftw3.h>
#endif

int main (int argc, char* argv[]) 
{
    int n1, n2, n, nk, i, j, k, nlags;
    float *data, wt;
    sf_complex **fft, *dataf; 
    sf_file inp, out;
#ifdef SF_HAS_FFTW
    fftwf_plan cfg=NULL, icfg=NULL;
#else
    kiss_fftr_cfg cfg=NULL, icfg=NULL;
#endif

    sf_init(argc, argv);

    inp = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");

    if (!sf_getint("nlags",&nlags)) nlags=100; /* number of lags */

    sf_putint(out,"n1",nlags);
    sf_putint(out,"n3",n2);

    nk = kiss_fft_next_fast_size((n1+1)/2)+1;
    n = 2*(nk-1);

    wt = 1.0/n;

    data = sf_floatalloc(n);
    dataf = sf_complexalloc(nk);
    fft = sf_complexalloc2(nk,n2);

#ifdef SF_HAS_FFTW
    cfg = fftwf_plan_dft_r2c_1d(n, data, (fftwf_complex *) dataf,
				FFTW_MEASURE);
    icfg = fftwf_plan_dft_c2r_1d(n, (fftwf_complex *) dataf, data,
				 FFTW_MEASURE);
    if (NULL == cfg || NULL == icfg) sf_error("FFT allocation failure");
#else
    cfg  = kiss_fftr_alloc(n,0,NULL,NULL);
    icfg  = kiss_fftr_alloc(n,1,NULL,NULL);
#endif
    
    for (i=0; i < n2; i++) {
	sf_floatread(data,n1,inp);
	for (k=n1; k < n; k++) {
	    data[k] = 0.0f;
	}
#ifdef SF_HAS_FFTW
	fftwf_execute(cfg);
#else
	kiss_fftr (cfg,data,(kiss_fft_cpx *) dataf);
#endif
	for (k=0; k < nk; k++) {
	    fft[i][k] = dataf[k];
	}
    }

/*************************************************
*                                                *
*  cross-correlate every trace with every other  *
*                                                *
*************************************************/

    for (i=0; i < n2; i++) {
	for (j=0; j < n2; j++) {
	    for (k=0; k < nk; k++) {
		dataf[k] = fft[i][k] * conjf(fft[j][k]);
	    }

#ifdef SF_HAS_FFTW
	    fftwf_execute(icfg);
#else
	    kiss_fftri(icfg,(kiss_fft_cpx *) dataf,data);
#endif

	    for (k=0; k < nlags; k++) {
		data[k] *= wt;
	    }

	    sf_floatwrite(data,nlags,out);
	}
    }

    exit(0);
}


