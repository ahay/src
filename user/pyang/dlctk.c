/* DLCT-K linear operator
discrete linear chirp transfrom for t-->f, fourier transform for x-->kx
*/
/*
  Copyright (C) 2013  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

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
#include <fftw3.h>

#include "dlct.h"
#include "dlctk.h"

fftwf_complex *q;
fftwf_plan fft,ifft;

static int n1, n2, L;
static float C, **sig;
static sf_complex **Sc;

void dlctk_init(int n1_, int n2_, int L_, float C_)
/*< allocate variables and initilize >*/
{
    n1=n1_;
    n2=n2_;
    L=L_;
    C=C_;
	
    q=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*n2);
    fft=fftwf_plan_dft_1d(n2, q, q, FFTW_FORWARD, FFTW_MEASURE);
    ifft=fftwf_plan_dft_1d(n2, q, q, FFTW_BACKWARD, FFTW_MEASURE);	
    sig=sf_floatalloc2(n1,n2);
    Sc=sf_complexalloc2(n1*L,n2);
}

void dlctk_lop(bool adj, bool add, int nm, int nd, sf_complex *mm, float *dd)
/*< DLCT-k linear operator >*/
{
	int i1, i2, il;

    	if (nm!=n1*n2*L) sf_error("%s: size mismatch: %d != %d",__FILE__, nm, n1*n2);
    	if (nd!=n1*n2) sf_error("%s: size mismatch: %d != %d",__FILE__, nd, n1*n2*L);

	if(adj){
		memcpy(sig[0],dd,nd*sizeof(float));
		if(!add) memset(mm,0,nm*sizeof(sf_complex));
	}else{
		memcpy(Sc[0],mm,nm*sizeof(sf_complex));
		if(!add) memset(dd,0,nd*sizeof(float));
	}

	if(adj){
	    	for(i2=0; i2<n2; i2++) forward_dlct(n1, L, C, sig[i2], Sc[i2]);
		for(i1=0; i1<n1; i1++) 
		for(il=0; il<L; il++)
		{
			for(i2=0; i2<n2; i2++) q[i2]=Sc[i2][i1+n1*il];
			fftwf_execute(fft);			
			for(i2=0; i2<n2; i2++) 
			{
				Sc[i2][i1+n1*il]=q[i2]/sqrtf(n2);
				mm[i2+n2*(i1+n1*il)]+=Sc[i2][i1+n1*il];
			}
		}
	}else{
		for(i1=0; i1<n1; i1++) 
		for(il=0; il<L; il++)
		{
			for(i2=0; i2<n2; i2++) q[i2]=Sc[i2][i1+n1*il];
			fftwf_execute(ifft);			
			for(i2=0; i2<n2; i2++) Sc[i2][i1+n1*il]=q[i2]/sqrtf(n2);
		}
	    	for(i2=0; i2<n2; i2++) inverse_dlct(n1, L, C, sig[i2], Sc[i2]);
		for(i2=0; i2<n2; i2++)
		for(i1=0; i1<n1; i1++) 
			dd[i1+n1*i2]+=sig[i2][i1];
	}
}

void dlctk_close()
/*< free the allovated variables >*/
{
    free(*sig); free(sig);
    free(*Sc); free(Sc);

    fftwf_free(q);
    fftwf_destroy_plan(fft);
    fftwf_destroy_plan(ifft);
}

