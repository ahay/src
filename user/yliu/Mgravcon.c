/* Continuation for gravity data by using FFT or intergral iteration */
/*
  Copyright (C) 2017 Jilin University
  
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

int main(int argc,char* argv[])
{
    bool verb;                 /* error check flag */
    bool itersymb;             /* iteration method or not */
    float a, z, b;             /* calculate parameter */  
    int n1, n2, nkx, nky, nw;  /* dimension info */
    float d1, o1, d2, o2, wkx, wky, dw;
    float shift;
    int i1, i2;  
    int niter, iter;           /* iteration number for iteration method */
    
    /* float **amp=NULL;       /\* amplititude spectrum *\/ */
    float **cof=NULL;          /* fft continuation filter */
    float **itf=NULL;          /* iteration continuation factor */
    float *dd=NULL, *dd2=NULL;           /* float data matrix   */
    kiss_fft_cpx **ff=NULL, **cc=NULL, ce;   /* complex data matrix */
    
    kiss_fft_cpx **ittmp=NULL;  /* iteration temp */  
    float *d;                          /* temp */
    kiss_fft_cpx *ctrace1, *ctrace2;   /* temp */
    kiss_fftr_cfg cfg1, icfg1;         
    kiss_fft_cfg cfg2, icfg2;
    
    sf_file in, out;
    
    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");
    
    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    if (!sf_getbool("iter",&itersymb)) itersymb = false;
    /* if y, perform iteration method */
  
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");

    if (!sf_histfloat(in,"d1",&d1)) d1 = 1. ;
    if (!sf_histfloat(in,"o1",&o1)) o1 = 0. ;

    if (!sf_histfloat(in,"d2",&d2)) d2 = 1. ;
    if (!sf_histfloat(in,"o2",&o2)) o2 = 0. ;

    /* determine wavenumber sampling (for real to complex) */
    nkx = 2*kiss_fft_next_fast_size((n1+1)/2); 
    nw = nkx/2+1;
    dw=1./(nkx*d1);

    /* detemine wavenumber sampling (for complex) */
    nky = kiss_fft_next_fast_size(n2*2);

    if(verb) sf_warning("wavenumber-determine done,kx=%d,ky=%d",nw,nky);  
    cfg1  = kiss_fftr_alloc(nkx,0,NULL,NULL);
    icfg1 = kiss_fftr_alloc(nkx,1,NULL,NULL);
    cfg2  = kiss_fft_alloc(nky,0,NULL,NULL);
    icfg2 = kiss_fft_alloc(nky,1,NULL,NULL);

    /* continuation deep */
    if (!sf_getfloat("z",&z)) sf_error("Need z= input");
  
    /* for iteration method */
    if (!sf_getint("niter",&niter)) niter = 0;

    /* continuation factor allocate memory */
    cof = sf_floatalloc2(nw,nky);
    itf = sf_floatalloc2(nw,nky);
    ittmp = (kiss_fft_cpx**)sf_complexalloc2(nw,nky);
  
    /* allocate memory for data compute */
    dd = sf_floatalloc(n1*n2);
    dd2 = sf_floatalloc(n1*n2);
    ff = (kiss_fft_cpx**)sf_complexalloc2(nw,nky);
    cc = (kiss_fft_cpx**)sf_complexalloc2(nw,nky);
  
    d = sf_floatalloc(nkx);
    ctrace1 = (kiss_fft_cpx*) sf_complexalloc(nw);    /* FFT1 trace */
    ctrace2 = (kiss_fft_cpx*) sf_complexalloc(nky);   /* FFT2 trace */
    if(verb) sf_warning("memory-allocate done");
  
    /* scaling for ifft */
    wkx = 1.0/nkx;
    wky = 1.0/nky;
  
    /* read in data  */
    sf_floatread(dd,n1*n2,in); 
  
    /* 2D-FFT start  */
    /* fft1 x->kx */
    for (i2=0; i2<n2; i2++) {
	for (i1=0; i1<nkx; i1++) {
	    if (i1<n1) {
		d[i1] = dd[i2*n1+i1];
	    } else {
		d[i1] = 0. ;                   /* padding with zero */
	    }
	}
	kiss_fftr(cfg1,d,ctrace1);
	for (i1=0; i1<nw; i1++)	{
	    ff[i2][i1] = i2%2?
		sf_cneg(ctrace1[i1]) : ctrace1[i1];  /* centerting  */
	}
	if (0. != o1) {
	    for (i1=0; i1 < nw; i1++) {
		shift = -2.0*SF_PI*i1*dw*o1;
		ce.r = cosf(shift);
		ce.i = sinf(shift);
		ff[i2][i1]=sf_cmul(ff[i2][i1],ce);
	    }
	}
    }
    for (i2=n2; i2<nky; i2++) {             /* padding with zero */
	for (i1=0; i1<nw; i1++)	{
	    ff[i2][i1].r = 0. ;
	    ff[i2][i1].i = 0. ;
	}
    }
    /* fft2 y->ky */
    for (i1=0; i1<nw; i1++) {
	kiss_fft_stride(cfg2,ff[0]+i1,ctrace2,nw);
	for (i2=0; i2<nky; i2++) {             /* transpose */
	    cc[i2][i1] = ctrace2[i2];           
	}
    }
    if(verb) sf_warning("fft2 done");
    /* 2D-FFT finish */
  
    /* FFT continuation filter/factor */
    for (i1=0; i1<nw; i1++) {
	for (i2=0; i2<nky; i2++) {
	    if (i2<nky/2) {
		b = sqrt((i1+1)*(i1+1)*0.01+(nky/2-i2)*(nky/2-i2)*0.01);
		a = exp(-1*z*b);
		cof[i2][i1] = a;
	    } else {
		b = sqrt((i1+1)*(i1+1)*0.01+abs(nky/2-i2)*abs(nky/2-i2)*0.01);
		a = exp(-1*z*b);
		cof[i2][i1] = a;
	    }	 
	}
    }
  
    if (itersymb) {     
	if(verb) sf_warning("Iterative continuation.");
	if (!niter) sf_error("Need niter= input");
	for (i1=0; i1<nw; i1++) {
	    for (i2=0; i2<nky; i2++) {
		itf[i2][i1] = 1.0-cof[i2][i1];  /* iteration factor */
		ittmp[i2][i1] = cc[i2][i1];     /* iteration initialization */
	    }
	}          
	for (iter=0; iter<niter; iter++) {
	    if(verb) sf_warning("iter %d of %d;",iter+1,niter);
	    for (i1=0; i1<nw; i1++) {
		for (i2=0; i1<nkx; i1++) {
		    cc[i2][i1] =
			sf_cadd(sf_crmul(cc[i2][i1],itf[i2][i1]),ittmp[i2][i1]);
		}
	    }
	}
	if(verb) sf_warning(".");
    } else { /* FFT continuation */
	if(verb) sf_warning("FFT continuation.");
	for (i1=0; i1<nw; i1++) {
	    for (i2=0; i2<nky; i2++) {
		cc[i2][i1] = sf_crmul(cc[i2][i1],cof[i2][i1]);
	    }
	}    
    }
  
    /* 2D-IFFT start */
    /* ifft1 ky->y */
    for (i1=0; i1<nw; i1++) {
	kiss_fft_stride(icfg2,cc[0]+i1,ctrace2,nw);
	for (i2=0; i2<nky; i2++) {
	    ff[i2][i1] = sf_crmul(ctrace2[i2],i2%2? -wky : wky);
	}
    }
    /* ifft2 kx->x */
    for (i2=0; i2<n2; i2++) {
	if (0. != o1) {
	    for (i1=0; i1 < nw; i1++) {
		shift = +2.0*SF_PI*i1*dw*o1;
		ce.r = cosf(shift);
		ce.i = sinf(shift);
		ff[i2][i1]=sf_cmul(ff[i2][i1],ce);
	    }
	}
	kiss_fftri(icfg1,ff[i2],d);
	for (i1=0; i1<n1; i1++) {
	    dd2[i2*n1+i1] = d[i1]*wkx;
	}
    }

    sf_floatwrite(dd2,n1*n2,out);

    free(*ff);   free(ff);   free(*cc); free(cc);
    free(*cof);  free(cof);  free(*itf);free(itf);
    free(*ittmp);free(ittmp);
    free(dd);    free(dd2);

    exit(0);  
}
