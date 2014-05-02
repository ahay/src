/* n-D Two-step POCS interpolation using a general Lp-norm optimization
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

  Reference: On analysis-based two-step interpolation methods for randomly 
	sampled seismic data, P Yang, J Gao, W Chen, Computers & Geosciences
*/
#include <rsf.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "pthresh.h"

int main(int argc, char* argv[])
{
    	bool verb;
    	int i, i1, i2, index, n1, n2, num, dim, n[SF_MAX_DIM], nw, iter, niter, nthr;
    	float thr, p, pclip, m;
    	//float t0, t1, beta;	
    	float *tdat, *thresh, *mask;
    	char *mode, key[7];
    	fftwf_complex *wdat, *wdat1;
    	fftwf_plan fft1, ifft1, fft2, ifft2;/* execute plan for FFT and IFFT */
    	sf_file in, out, Fmask;/* mask and I/O files*/ 


    	sf_init(argc,argv);	/* Madagascar initialization */
#ifdef _OPENMP
    	omp_init(); 	/* initialize OpenMP support */
#endif

    	/* setup I/O files */
    	in=sf_input("in");	/* read the data to be interpolated */
    	out=sf_output("out"); 	/* output the reconstructed data */
    	Fmask=sf_input("mask");  	/* read the (n-1)-D mask for n-D data */
 
    	if(!sf_getbool("verb",&verb))    	verb=false;
    	/* verbosity */
    	if (!sf_getint("niter",&niter)) 	niter=100;
    	/* total number iterations */
    	if (!sf_getfloat("pclip",&pclip)) 	pclip=10.;
    	/* starting data clip percentile (default is 99)*/
    	if (pclip <=0. || pclip > 100.)
	sf_error("pclip=%g should be > 0 and <= 100",pclip);
    	if ( !(mode=sf_getstring("mode")) ) mode = "exp";
    	/* thresholding mode: 'hard', 'soft','pthresh','exp';
	'hard', hard thresholding;	'soft', soft thresholding; 
	'pthresh', generalized quasi-p; 'exp', exponential shrinkage */
    	if (!sf_getfloat("p",&p)) 		p=0.35;
    	/* norm=p, where 0<p<=1 */;
    	if (strcmp(mode,"soft") == 0) 		p=1;
    	else if (strcmp(mode,"hard") == 0) 	p=0;

    	/* dimensions */
   	for (i=0; i < SF_MAX_DIM; i++) {
		snprintf(key,3,"n%d",i+1);
		if (!sf_getint(key,n+i) && 
		    (NULL == in || !sf_histint(in,key,n+i))) break;
		/*( n# size of #-th axis )*/  
		sf_putint(out,key,n[i]);
    	}

    	if (0==i) sf_error("Need n1=");
    	dim=i;

    	n1=n[0];
    	n2=sf_leftsize(in,1);
	nw=n1/2+1;
	num=nw*n2;//total number of elements in frequency domain
 
    	/* allocate data and mask arrays */
	thresh=(float*)malloc(nw*n2*sizeof(float));
    	tdat=(float*)fftwf_malloc(n1*n2*sizeof(float)); // data in time domain
    	wdat=(fftwf_complex*)fftwf_malloc(nw*n2*sizeof(fftwf_complex));// data in frequency domain
    	wdat1=(fftwf_complex*)fftwf_malloc(nw*n2*sizeof(fftwf_complex));// data in frequency domain
    	fft1=fftwf_plan_many_dft_r2c(1, &n1, n2, tdat, &n1, 1, n1, wdat, &n1, 1, nw, FFTW_MEASURE);	
   	ifft1=fftwf_plan_many_dft_c2r(1, &n1, n2, wdat, &n1, 1, nw, tdat, &n1, 1, n1, FFTW_MEASURE);
	fft2=fftwf_plan_many_dft(dim-1, &n[1], nw, wdat1, &n[1], nw, 1, wdat1, &n[1], nw, 1, FFTW_FORWARD,FFTW_MEASURE);
	ifft2=fftwf_plan_many_dft(dim-1, &n[1], nw, wdat1, &n[1], nw, 1, wdat1, &n[1], nw, 1, FFTW_BACKWARD,FFTW_MEASURE);

	/* initialization */
    	sf_floatread(tdat, n1*n2, in);
    	if (NULL != sf_getstring("mask")){
		mask=sf_floatalloc(n2);
		sf_floatread(mask,n2,Fmask);
    	}

	// transform the data from time domain to frequency domain: tdat-->wdat 
	fftwf_execute(fft1);
	for(i=0; i<num; i++) wdat[i]/=sqrtf(n1);
	memset(wdat1,0,num*sizeof(fftwf_complex));

    	for(iter=1; iter<=niter; iter++)
    	{
		fftwf_execute(fft2);
	#ifdef _OPENMP
	#pragma omp parallel for default(none)	\
		private(i)			\
		shared(wdat1,num,n2)
	#endif
		for(i=0; i<num; i++) wdat1[i]/=sqrtf(n2);

		// perform hard thresholding
	#ifdef _OPENMP
	#pragma omp parallel for default(none)	\
		private(i)			\
		shared(thresh,wdat1,num)
	#endif
		for(i=0; i<num; i++)	thresh[i]=cabsf(wdat1[i]);

	   	nthr = 0.5+num*(1.-0.01*pclip); 
	    	if (nthr < 0) nthr=0;
	    	if (nthr >= num) nthr=num-1;
		thr=sf_quantile(nthr,num,thresh);
		//thr*=powf(0.01,(iter-1.0)/(niter-1.0));

	#ifdef _OPENMP
	#pragma omp parallel for default(none)	\
		private(i)			\
		shared(wdat1,num,thr)
	#endif
		for(i=0; i<num; i++) wdat1[i]*=(cabsf(wdat1[i])>thr?1.:0.);

		fftwf_execute(ifft2);
	#ifdef _OPENMP
	#pragma omp parallel for default(none)	\
		private(i)			\
		shared(wdat1,num,n2)
	#endif
		for(i=0; i<num; i++) wdat1[i]/=sqrtf(n2);		

		// d_rec = d_obs+(1-M)*A T{ At(d_rec) } 
	#ifdef _OPENMP
	#pragma omp parallel for default(none)	\
		private(i1,i2,index,m)		\
		shared(mask,wdat,wdat1,nw,n2)
	#endif
		for(i2=0; i2<n2; i2++)
		for(i1=0; i1<nw; i1++)
		{ 
			m=(mask[i2])?1:0;
			index=i1+nw*i2;
			wdat1[index]=wdat[index]+(1.-m)*wdat1[index];
		}

		if (verb)    sf_warning("iteration %d;",iter);
    	}

	// transform the data from frequency domain to time domain: wdat-->tdat
	memcpy(wdat, wdat1, num*sizeof(fftwf_complex));
	fftwf_execute(ifft1);
	for(i=0; i<n1*n2; i++) tdat[i]/=sqrtf(n1);
    	sf_floatwrite(tdat, n1*n2, out);/* output reconstructed seismograms */

	free(thresh);
	fftwf_free(tdat);
	fftwf_free(wdat);
	fftwf_free(wdat1);	

    	exit(0);
}
