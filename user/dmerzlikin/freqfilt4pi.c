/* Frequency-domain filtering in 2-D */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include <math.h>

#include <rsf.h>
/*^*/

#include "freqfilt4pi.h"
#include "Faddeeva.h"

static int nfft, nw, m1, m2;
static kiss_fft_cpx *ctrace, *ctrace2, **fft;
static float *trace/*, **shape */;
static sf_complex **shape;
kiss_fftr_cfg tfor, tinv;
kiss_fft_cfg  xfor, xinv;

void freqfilt4pi_init(int n1, int n2 /* data dimensions */, 
		    int nw1        /* number of frequencies */)
/*< initialize >*/
{
    m1 = n1;
    nw = nw1;
    m2 = n2;
    nfft = 2*kiss_fft_next_fast_size((n1+1)/2);

    tfor = kiss_fftr_alloc(nfft,0,NULL,NULL);
    tinv = kiss_fftr_alloc(nfft,1,NULL,NULL);
    xfor = kiss_fft_alloc(n2,0,NULL,NULL);
    xinv = kiss_fft_alloc(n2,1,NULL,NULL);
    if (NULL == tfor || NULL == tinv || NULL == xfor || NULL == xinv)
	sf_error("%s: KISS FFT allocation error",__FILE__);

    trace = sf_floatalloc(nfft);
    ctrace = (kiss_fft_cpx*) sf_complexalloc(nw);
    ctrace2 = (kiss_fft_cpx*) sf_complexalloc(n2);
    fft = (kiss_fft_cpx**) sf_complexalloc2(nw,n2);
}

void freqfilt4pi_set(/*float*/ sf_complex **filt)
/*< set the filter >*/
{
    shape = filt;
}

void freqfilt4pi_close(void) 
/*< free allocated storage >*/
{
    free (tfor);
    free (tinv);
    free (xfor);
    free (xinv);
    free (trace);
    free (ctrace);
    free (ctrace2);
    free (fft[0]);
    free (fft);
}

void freqfilt4pi_spec (const float* x /* input */, float** y /* spectrum */) 
/*< compute 2-D spectrum >*/
{
    int ik, iw;

    for (ik=0; ik < m2; ik++) {
	for (iw=0; iw < m1; iw++) {
	    trace[iw] = x[ik*m1+iw];
	}	
	for (iw=m1; iw < nfft; iw++) {
	    trace[iw]=0.;
	}

	kiss_fftr (tfor,trace,ctrace);
	for (iw=0; iw < nw; iw++) {
	    fft[ik][iw] = ik%2? sf_cneg(ctrace[iw]): ctrace[iw];
	}
    }

    for (iw=0; iw < nw; iw++) {
	kiss_fft_stride(xfor,fft[0]+iw,ctrace2,nw);
	for (ik=0; ik < m2; ik++) {
	    y[iw][ik] = sf_cabsf(ctrace2[ik]); /* transpose */
	}
    }
}

void freqfilt4pi_lop (bool adj, bool add, int nx, int ny, float* x, float* y) 
/*< linear filtering operator >*/
{
    int iw, ik;
    kiss_fft_cpx temp;

    sf_adjnull(adj,add,nx,ny,x,y);

    for (ik=0; ik < m2; ik++) {
	for (iw=0; iw < m1; iw++) {
	    trace[iw] = adj? y[ik*m1+iw]: x[ik*m1+iw];
	
		/*if(trace[iw] == 0.0){
		
			sf_warning("trace[%d] = %f",iw,trace[iw]);
			
		}*/
	
	}
	for (iw=m1; iw < nfft; iw++) {
	    trace[iw]=0.;
	}

	kiss_fftr (tfor,trace,ctrace);
	for (iw=0; iw < nw; iw++) {
	    fft[ik][iw] = ik%2? sf_cneg(ctrace[iw]): ctrace[iw];
	}
    }

    for (iw=0; iw < nw; iw++) {
	kiss_fft_stride(xfor,fft[0]+iw,ctrace2,nw);

	for (ik=0; ik < m2; ik++) {
	
		//double creal( double complex z );
	
	    //transform to kiss fft cpx
	    //creal are double complex functions - what should we
	    //do when double complex is not supported???  	
		if (adj){
			temp.r = creal(shape[iw][ik]);
			temp.i = (-1.0)*cimag(shape[iw][ik]);
		} else {
			temp.r = creal(shape[iw][ik]);
			temp.i = cimag(shape[iw][ik]);
	    }
	    
	    ctrace2[ik] = sf_cmul(ctrace2[ik],temp);
	
	}

	kiss_fft(xinv,ctrace2,ctrace2);

	for (ik=0; ik < m2; ik++) {
	    fft[ik][iw] = ik%2? sf_cneg(ctrace2[ik]): ctrace2[ik];
	}
    }

    for (ik=0; ik < m2; ik++) {
	kiss_fftri (tinv,fft[ik],trace);

	for (iw=0; iw < m1; iw++) {	  
	    if (adj) {
		x[ik*m1+iw] += trace[iw];
	    } else {
		y[ik*m1+iw] += trace[iw];
	    }
	}
    }
}

/* 	$Id: freqfilt4pi.c 2262 2006-09-15 04:50:52Z sfomel $	 */
