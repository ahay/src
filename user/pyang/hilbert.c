/* Hilbert transform of seismic traces 
Note that the 1st dimension and 2nd dimension have been flipped. 2nd: time
*/
/*
  Copyright (C) 2014  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

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
#include <complex.h>

#include "hilbert.h"

#ifdef SF_HAS_FFTW
#include <fftw3.h>

fftwf_plan fftn, ifftn;/* execute plan for FFT and IFFT */
fftwf_complex *fftout;
float *fftin;
int n1g, n2t;

void hilbert_init(int ng, int nt)
/*< initialize hilbert transform: n1=ng; n2=nt; >*/
{
	n1g=ng;
	n2t=nt;
	fftin=(float*)fftwf_malloc(ng*nt*sizeof(float));
	fftout=(fftwf_complex*)fftwf_malloc(ng*nt*sizeof(fftwf_complex));
	fftn=fftwf_plan_many_dft_r2c(1, &nt, ng, fftin, &nt, ng, 1,fftout, &nt, ng, 1,FFTW_MEASURE);
	ifftn=fftwf_plan_many_dft(1, &nt, ng, fftout, &nt, ng, 1,fftout, &nt, ng, 1,FFTW_BACKWARD,FFTW_MEASURE);
}


void hilbert_trans(float *in, sf_complex *out)
/*< hilbert transform: adobs is an analytical signal >*/
{
    	int ig, it;

	memcpy(fftin, in, n1g*n2t*sizeof(float));
	fftwf_execute(fftn);/* FFT */
	for(ig=0; ig<n1g; ig++)
	{
		/* positive freqs: x 2 */
		for(it=1; it<n2t/2; it++) 	fftout[ig+it*n1g]*=2.;
		/* negtive freqs: 0 */
		for(it=n2t/2+1; it<n2t; it++)	fftout[ig+it*n1g]=0;
	}
	fftwf_execute(ifftn);/* IFFT */
	for(it=0; it<n2t; it++)
	for(ig=0; ig<n1g; ig++)
		fftout[ig+it*n1g]/=n1g*n2t;
	memcpy(out, fftout, n1g*n2t*sizeof(sf_complex));
}

void hilbert_close()
/*< free variables >*/
{
    fftwf_free(fftin);
    fftwf_free(fftout);
    fftwf_destroy_plan(fftn);
    fftwf_destroy_plan(ifftn);
}
#endif
