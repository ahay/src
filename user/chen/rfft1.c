/* 1D real data fft interface */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
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

typedef struct tag_rfft1
{
    int nt;		/*	real data */
    int nfft;	        /*	fft size with zero padding */
    int nk;		/*	number of frequencies */
    kiss_fftr_cfg fwd;
    kiss_fftr_cfg inv;
    float 		*inbuf;
    float 		*outbuf;
}rfft1;

void *sf_rfft1_init(int *nt, int *nk, int pad)
/*< initialize: return a handle to avoid public data >*/
{
    rfft1 *p;
    int n1;

    p = (rfft1*) malloc(sizeof(rfft1));
    if(p==NULL) return NULL;

    n1 = kiss_fft_next_fast_size( (*nt+1)/2 );
    p->nt = *nt;
    if(pad) p->nfft = n1*4;
    else	p->nfft = n1*2;
    p->nk = p->nfft/2+1;

    p->inbuf = sf_floatalloc(p->nk*4);
    if(p->inbuf==NULL) {free(p);  return NULL;}
    p->outbuf = p->inbuf + p->nk*2;

    p->fwd = kiss_fftr_alloc(p->nfft, 0, NULL, NULL);
    p->inv = kiss_fftr_alloc(p->nfft, 1, NULL, NULL);

    *nt = p->nfft;
    *nk = p->nk;
    return p;
}

void sf_rfft1(void *h, float *in, sf_complex *out)
/*< forward fft: only nt accesed >*/
{
    int i1;
    rfft1 *p;

    p = (rfft1*) h;
    for(i1=0; i1<p->nt; i1++) p->inbuf[i1] = in[i1];
    for(; i1<p->nk*2; i1++) p->inbuf[i1] = 0.0;

    kiss_fftr (p->fwd, p->inbuf, (kiss_fft_cpx*) p->outbuf);

    for(i1=0; i1<p->nk; i1++) 
	out[i1] = ((sf_complex*)p->outbuf)[i1];
}


void sf_rifft1(void *h, sf_complex *in, float *out)
/*< inverse fft >*/
{
    int i1;
    rfft1 *p;

    p = (rfft1*) h;
    for(i1=0; i1<p->nk; i1++) 
	((sf_complex*)p->inbuf)[i1] = in[i1];

    kiss_fftri (p->inv, (kiss_fft_cpx*) p->inbuf, p->outbuf);

    /* only nt data returned, so we do not need alloc memory in the parent function */
    for(i1=0; i1<p->nt; i1++) 
	out[i1] = p->outbuf[i1]/p->nfft;
}

void sf_rfft1_close(void *h)
/*< release memory >*/
{
    kiss_fftr_free(((rfft1*)h)->fwd);
    kiss_fftr_free(((rfft1*)h)->inv);
    free((rfft1*)h);
}


