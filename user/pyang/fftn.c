/* n-d fft operator 
Note: The adjoint is made as the same as inverse by normalization!
*/

/*
  Copyright (C) 2013 Pengliang Yang, Xi'an Jiaotong University, UT Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.h>

#ifdef SF_HAS_FFTW
#include <fftw3.h>

#include "fftn.h"

static int num;
fftwf_plan fftn, ifftn;/* execute plan for FFT and IFFT */
fftwf_complex *tmp;

void fftn_init(int rank, int *n)
/*< initialize >*/
{
	int i;
	num=1;
	for(i=0; i<rank; i++) num*=n[i];
    	tmp=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*num);
    	fftn=fftwf_plan_dft(rank, n, tmp, tmp, FFTW_FORWARD, FFTW_MEASURE);	
   	ifftn=fftwf_plan_dft(rank, n, tmp, tmp, FFTW_BACKWARD, FFTW_MEASURE);
}

void fftn_lop (bool adj, bool add, int nx, int ny, sf_complex *xx, sf_complex *yy)
/*< linear operator >*/
{
    int i;

    if (num!=nx) sf_error("%s: size mismatch: %d != %d",__FILE__, num, nx);
    if (num!=ny) sf_error("%s: size mismatch: %d != %d",__FILE__, num, ny);

    sf_cadjnull (adj, add, nx, ny, xx, yy);
    for(i=0; i<num; i++) tmp[i] = adj? yy[i]: xx[i];
  
    if (adj){
 	fftwf_execute(fftn);	
    	for(i=0; i<num; i++) xx[i]+=tmp[i]/sqrtf(num);
    }else{
	fftwf_execute(ifftn);
    	for(i=0; i<num; i++) yy[i]+=tmp[i]/sqrtf(num);
    }
}

void fftn_close(void)
/*< free allocated variables>*/
{
    fftwf_free(tmp);
    fftwf_destroy_plan(fftn);
    fftwf_destroy_plan(ifftn);
}


#endif
