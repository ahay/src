/* 3-D FFT interface */
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

static bool cmplx;
static int n1, n2, n3, nk;
static float wt;
static kiss_fftr_cfg cfg, icfg;
static kiss_fft_cfg cfg1, icfg1, cfg2, icfg2, cfg3, icfg3;
static kiss_fft_cpx ***tmp, *ctrace1, *ctrace2, *ctrace3;
static sf_complex *trace1, *trace2, *trace3;

int fft3_init(bool cmplx1        /* if complex transform */,
	      int pad1           /* padding on the first axis */,
	      int nx,   int ny,   int nz   /* input data size */, 
	      int *nx2, int *ny2, int *nz2 /* padded data size */)
/*< initialize >*/
{
    int i2, i3;

    cmplx = cmplx1;

    /* axis 1 */

    if (cmplx) {
	nk = n1 = kiss_fft_next_fast_size(nx*pad1);
	
	cfg1  = kiss_fft_alloc(n1,0,NULL,NULL);
	icfg1 = kiss_fft_alloc(n1,1,NULL,NULL);

	trace1 = sf_complexalloc(n1);
	ctrace1 = (kiss_fft_cpx *) trace1;
    } else {
	nk = kiss_fft_next_fast_size(pad1*(nx+1)/2)+1;
	n1 = 2*(nk-1);

	cfg  = kiss_fftr_alloc(n1,0,NULL,NULL);
	icfg = kiss_fftr_alloc(n1,1,NULL,NULL);
    }

    *nx2 = n1;

    /* axis 2 */

    n2 = kiss_fft_next_fast_size(ny);

    cfg2  = kiss_fft_alloc(n2,0,NULL,NULL);
    icfg2 = kiss_fft_alloc(n2,1,NULL,NULL);

    *ny2 = n2;

    trace2 = sf_complexalloc(n2);
    ctrace2 = (kiss_fft_cpx *) trace2;

    /* axis 3 */

    n3 = kiss_fft_next_fast_size(nz);

    cfg3  = kiss_fft_alloc(n3,0,NULL,NULL);
    icfg3 = kiss_fft_alloc(n3,1,NULL,NULL);

    *nz2 = n3;

    trace3 = sf_complexalloc(n3);
    ctrace3 = (kiss_fft_cpx *) trace3;

    /* --- */

    tmp = (kiss_fft_cpx***) sf_alloc (n3,sizeof(kiss_fft_cpx**));
    tmp[0] = (kiss_fft_cpx**) sf_alloc (n2*n3,sizeof(kiss_fft_cpx*));
    tmp[0][0] = (kiss_fft_cpx*) sf_alloc (nk*n2*n3,sizeof(kiss_fft_cpx));

    for (i2=1; i2 < n2*n3; i2++) {
	tmp[0][i2] = tmp[0][0]+i2*nk;
    }

    for (i3=1; i3 < n3; i3++) {
	tmp[i3] = tmp[0]+i3*n2;
    }


    wt =  1.0/(n3*n2*n1);

    return (nk*n2*n3);
}

void fft3(float ***inp      /* [n1*n2] */, 
	  sf_complex ***out /* [nk*n2] */)
/*< 3-D FFT >*/
{
    int i1, i2, i3;

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    if (cmplx) {
		for (i1=0; i1 < n1; i1++) {
		    trace1[i1] = sf_cmplx(inp[i3][i2][i1],0.);
		    if (i1%2) ctrace1[i1] = sf_cneg(ctrace1[i1]);
		}
		kiss_fft_stride(cfg1,ctrace1,tmp[i3][i2],1);
	    } else {
		kiss_fftr (cfg,inp[i3][i2],tmp[i3][i2]);
	    }
	}
    }

    /* FFT over second axis */
    for (i3=0; i3 < n3; i3++) {
	for (i2=1; i2<n2; i2+=2) {
	    for (i1=0; i1<nk; i1++) {
		tmp[i3][i2][i1] = sf_cneg(tmp[i3][i2][i1]);
	    }
	}

    
	for (i1=0; i1 < nk; i1++) {
	    kiss_fft_stride(cfg2,tmp[i3][0]+i1,ctrace2,nk);
	    for (i2=0; i2 < n2; i2++) {
		tmp[i3][i2][i1]=ctrace2[i2];
	    }
	}
    }

    /* FFT over third axis */

    for (i3=1; i3 < n3; i3+=2) {
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1<nk; i1++) {
		tmp[i3][i2][i1] = sf_cneg(tmp[i3][i2][i1]);
	    }
	}
    }

    
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < nk; i1++) {
	    kiss_fft_stride(cfg3,tmp[0][0]+i2*nk+i1,ctrace3,nk*n2);
		
	    /* Transpose */
	    for (i3=0; i3<n3; i3++) {
		out[i3][i2][i1] = trace3[i3];
	    }
	}
    }    
}


void ifft3(float ***out      /* [n1*n2] */, 
	   sf_complex ***inp /* [nk*n2] */)
/*< 3-D inverse FFT >*/
{
    int i1, i2, i3;

    /* IFFT over third axis */

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < nk; i1++) {
	    kiss_fft_stride(icfg3,(kiss_fft_cpx *) (inp[0][0]+i2*nk+i1),ctrace3,nk*n2);
	    
	    for (i3=0; i3<n3; i3++) {
		tmp[i3][i2][i1] = i3%2? sf_cneg(ctrace3[i3]):ctrace3[i3];
	    }
	}
    }
    
    /* IFFT over second axis */

    for (i3=0; i3 < n3; i3++) {
	for (i1=0; i1 < nk; i1++) {
	    kiss_fft_stride(icfg2,tmp[i3][0]+i1,ctrace2,nk);
		
	    for (i2=0; i2<n2; i2++) {
		tmp[i3][i2][i1] = sf_crmul(ctrace2[i2],i2%2? -wt: wt);
	    }
	}
    }

    /* IFFT over first axis */

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    if (cmplx) {
		kiss_fft_stride(icfg1,tmp[i3][i2],ctrace1,1);
		for (i1=0; i1 < n1; i1++) {
		    if (i1%2) ctrace1[i1] = sf_cneg(ctrace1[i1]);
		   // out[(i3*n2+i2)*n1+i1] = crealf(trace1[i1]);
                    out[i3][i2][i1] = crealf(trace1[i1]);
		}
	    } else {
		kiss_fftri(icfg,tmp[i3][i2],out[i3][i2]);
	    }
	}
    }
}
