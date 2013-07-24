/* 2D real data fft interface */

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

typedef struct tag_rfft2
{
    int nt1;		/*	real data */
    int nt2;		/*	real data */
    int nfft1;	        /*	fft size with zero padding */
    int nfft2;	        /*	fft size with zero padding */
    int nk1;		/*	number of frequencies */
    int nk2;		/*	number of frequencies */
    kiss_fftr_cfg fwd1;
    kiss_fftr_cfg inv1;
    kiss_fft_cfg fwd2;
    kiss_fft_cfg inv2;
    float 		*v1;	/* nfft1 */
    sf_complex	*v2;	/* nk1 */
    sf_complex	*v3;	/* nfft2 */
    sf_complex	**u1;	/* nk1*nfft2 */
}rfft2;

void *sf_rfft2_init(int *nt1, int *nt2, int *nk1, int pad)
/*< initialize: return a handle to avoid public data >*/
{
    rfft2 *p;
    int n1, n2;

    p = (rfft2*) malloc(sizeof(rfft2));
    if(p==NULL) return NULL;

    n1 = kiss_fft_next_fast_size( (nt1[0]+1)/2 );
    p->nt1 = nt1[0];

    n2 = kiss_fft_next_fast_size( (nt2[0]+1)/2 );
    p->nt2 = nt2[0];

    if(pad)
    {
	p->nfft1 = n1*4;
	p->nfft2 = n2*4;
    }else{
	p->nfft1 = n1*2;
	p->nfft2 = n2*2;
    }
    p->nk1 = p->nfft1/2+1;
    p->nk2 = p->nfft2;

    p->v1 = sf_floatalloc(p->nfft1);
    if(p->v1==NULL) {free(p);  return NULL;}
    p->v2 = sf_complexalloc(p->nk1);
    if(p->v2==NULL) {free(p); free(p->v1); return NULL;}
    p->v3 = sf_complexalloc(p->nk2);
    if(p->v3==NULL) {free(p); free(p->v1); free(p->v2); return NULL;}

    p->u1 = sf_complexalloc2(p->nk2, p->nk1);
    if(p->u1==NULL) 
    {
	free(p); 
	free(p->v1); 
	free(p->v2); 
	free(p->v3); 
	return NULL;
    }

    p->fwd1 = kiss_fftr_alloc(p->nfft1, 0, NULL, NULL);
    p->inv1 = kiss_fftr_alloc(p->nfft1, 1, NULL, NULL);

    p->fwd2 = kiss_fft_alloc(p->nfft2, 0, NULL, NULL);
    p->inv2 = kiss_fft_alloc(p->nfft2, 1, NULL, NULL);

    nt1[0] = p->nfft1;
    nt2[0] = p->nfft2;
    nk1[0] = p->nk1;

    return p;
}


void sf_rfft2(void *h, float **in, sf_complex **out)
/*< forward 2D fft: F(nt , nx) ==> C(nk1 , nfft2) >*/
{
    int i1, i2;
    rfft2 *p;

    p = (rfft2*) h;
    for(i2=0; i2<p->nt2; i2++)
    {
	for(i1=0; i1<p->nt1; i1++) p->v1[i1] = in[i2][i1];
	for(; i1<p->nfft1; i1++) p->v1[i1] = 0.0;
	kiss_fftr (p->fwd1, p->v1, (kiss_fft_cpx*) p->v2);
	for(i1=0; i1<p->nk1; i1++)	p->u1[i1][i2] = p->v2[i1];
    }
    for(i1=0; i1<p->nk1; i1++)
    {
	for(i2=p->nt2; i2<p->nfft2; i2++)
	    p->u1[i1][i2] = 0.0;
	kiss_fft (p->fwd2, (kiss_fft_cpx*) p->u1[i1], (kiss_fft_cpx*) p->v3);
	for(i2=0; i2<p->nfft2; i2++)	out[i2][i1] = p->v3[i2];
    }

}



void sf_rifft2(void *h, sf_complex **in, float **out)
/*< inverse 2D fft: C(nk1 , nfft2) ==> F(nt , nx) >*/
{
    int i1, i2;
    rfft2 *p;

    p = (rfft2*) h;
    for(i1=0; i1<p->nk1; i1++)
    {
	for(i2=0; i2<p->nfft2; i2++)	p->v3[i2] = in[i2][i1];
	kiss_fft(p->inv2, (kiss_fft_cpx*) p->v3, (kiss_fft_cpx*)p->u1[i1]);
    }
    for(i2=0; i2<p->nfft2; i2++)
    {
	for(i1=0; i1<p->nk1; i1++)	p->v2[i1] = p->u1[i1][i2];
	kiss_fftri (p->inv1, (kiss_fft_cpx*) p->v2, p->v1);
	if(i2 < p->nt2)
	    for(i1=0; i1 < p->nt1; i1++)
		out[i2][i1] = p->v1[i1]/(p->nfft1*p->nfft2);
    }

}


void sf_rfft2_close(void *h)
/*< release memory >*/
{
    kiss_fftr_free(((rfft2*)h)->fwd1);
    kiss_fftr_free(((rfft2*)h)->inv1);
    kiss_fft_free(((rfft2*)h)->fwd2);
    kiss_fft_free(((rfft2*)h)->inv2);
    free(((rfft2*)h)->v1); 
    free(((rfft2*)h)->v2); 
    free(((rfft2*)h)->v3); 
    free((rfft2*)h);
}


