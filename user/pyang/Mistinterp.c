/* n-D IST interpolation using a general Lp-norm optimization
Note: acquistion geometry illustrated by mask operator
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
    	float thr, pclip;
    	float *dobs_t, *thresh, *mask;
    	char key[7];
    	fftwf_complex *mm, *dd, *dobs;
    	fftwf_plan fft1, ifft1, fftrem, ifftrem;/* execute plan for FFT and IFFT */
    	sf_file in, out, Fmask;/* mask and I/O files*/ 


    	sf_init(argc,argv);/* Madagascar initialization */
#ifdef _OPENMP
    	omp_init(); 	/* initialize OpenMP support */
#endif

    	/* setup I/O files */
    	in=sf_input("in");	/* read the data to be interpolated */
    	out=sf_output("out"); 	/* output the reconstructed data */
    	Fmask=sf_input("mask");	/* read the (n-1)-D mask for n-D data */
 
    	if(!sf_getbool("verb",&verb))    	verb=false;
    	/* verbosity */
    	if (!sf_getint("niter",&niter)) 	niter=100;
    	/* total number iterations */
    	if (!sf_getfloat("pclip",&pclip)) 	pclip=10.;
    	/* starting data clip percentile (default is 99)*/
    	if (pclip <=0. || pclip > 100.)
	sf_error("pclip=%g should be > 0 and <= 100",pclip);

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
	num=nw*n2;/* data: total number of elements in frequency domain */
 
    	/* allocate data and mask arrays */
	thresh=(float*)malloc(nw*n2*sizeof(float));
    	dobs_t=(float*)fftwf_malloc(n1*n2*sizeof(float)); /* observed data in time domain*/
    	dobs=(fftwf_complex*)fftwf_malloc(nw*n2*sizeof(fftwf_complex)); /* observed data in frequency domain*/
    	dd=(fftwf_complex*)fftwf_malloc(nw*n2*sizeof(fftwf_complex));
    	mm=(fftwf_complex*)fftwf_malloc(nw*n2*sizeof(fftwf_complex));

	/* FFT for 1st dimension */
    	fft1=fftwf_plan_many_dft_r2c(1, &n1, n2, dobs_t, &n1, 1, n1, dobs, &n1, 1, nw, FFTW_MEASURE);	
   	ifft1=fftwf_plan_many_dft_c2r(1, &n1, n2, dobs, &n1, 1, nw, dobs_t, &n1, 1, n1, FFTW_MEASURE);
	/* FFT for remaining dimensions */
	fftrem=fftwf_plan_many_dft(dim-1, &n[1], nw, dd, &n[1], nw, 1, dd, &n[1], nw, 1, FFTW_FORWARD, FFTW_MEASURE);
	ifftrem=fftwf_plan_many_dft(dim-1, &n[1], nw, dd, &n[1], nw, 1, dd, &n[1], nw, 1, FFTW_BACKWARD, FFTW_MEASURE);

	/* initialization */
    	sf_floatread(dobs_t, n1*n2, in);
    	if (NULL != sf_getstring("mask")){
		mask=sf_floatalloc(n2);
		sf_floatread(mask,n2,Fmask);
    	}

	/* transform the data from time domain to frequency domain: dobs_t-->dobs */
	fftwf_execute(fft1);
	for(i=0; i<num; i++) dobs[i]/=sqrtf(n1);
	memset(mm,0,num*sizeof(fftwf_complex));

    	for(iter=1; iter<=niter; iter++)
    	{
		/* dd<-- A mm^k */
		memcpy(dd, mm, num*sizeof(fftwf_complex));
		fftwf_execute(ifftrem);
	#ifdef _OPENMP
	#pragma omp parallel for default(none)	\
		private(i)			\
		shared(dd,num,n2)
	#endif
		for(i=0; i<num; i++) dd[i]/=sqrtf(n2);

		/* apply mask: dd<-- dobs-M dd */
	#ifdef _OPENMP
	#pragma omp parallel for default(none)	\
		private(i1,i2,index)		\
		shared(mask,dobs,dd,nw,n2)
	#endif
		for(i2=0; i2<n2; i2++)
		for(i1=0; i1<nw; i1++)
		{ 
			index=i1+nw*i2;
			dd[index]=dobs[index]-mask[i2]*dd[index];
		}

		/* apply adjoint mask: dd<-- M^* dd (M^*=M) */
	#ifdef _OPENMP
	#pragma omp parallel for default(none)	\
		private(i1,i2,index)		\
		shared(mask,dd,nw,n2)
	#endif
		for(i2=0; i2<n2; i2++)
		for(i1=0; i1<nw; i1++)
		{ 
			index=i1+nw*i2;
			dd[index]=mask[i2]*dd[index];
		}

		/* mm^k += A^* dd */
		fftwf_execute(fftrem);
	#ifdef _OPENMP
	#pragma omp parallel for default(none)	\
		private(i)			\
		shared(dd,mm,num,n2)
	#endif
		for(i=0; i<num; i++) mm[i]+=dd[i]/sqrtf(n2);
		

		/* perform thresholding */
	#ifdef _OPENMP
	#pragma omp parallel for default(none)	\
		private(i)			\
		shared(thresh,mm,num)
	#endif
		for(i=0; i<num; i++)	thresh[i]=cabsf(mm[i]);


	   	nthr = 0.5+num*(1.-0.01*pclip); 
	    	if (nthr < 0) nthr=0;
	    	if (nthr >= num) nthr=num-1;
		thr=sf_quantile(nthr,num,thresh);
		/* thr*=powf(0.01,(iter-1.0)/(niter-1.0)); */

	#ifdef _OPENMP
	#pragma omp parallel for default(none)	\
		private(i)			\
		shared(mm,num,thr)
	#endif
		for(i=0; i<num; i++) mm[i]*=(cabsf(mm[i])>thr?1.:0.);

		if (verb)    sf_warning("iteration %d;",iter);
    	}

	/* transform the data from frequency domain to time domain: dobs-->dobs_t */
	memcpy(dd, mm, num*sizeof(fftwf_complex));
	fftwf_execute(ifftrem);
	for(i=0; i<num; i++) dd[i]/=sqrtf(n2);
	memcpy(dobs, dd, num*sizeof(fftwf_complex));
	fftwf_execute(ifft1);
	for(i=0; i<n1*n2; i++) dobs_t[i]/=sqrtf(n1);
    	sf_floatwrite(dobs_t, n1*n2, out);/* output reconstructed seismograms */

	free(thresh);
	fftwf_free(dobs_t);
	fftwf_free(dobs);
	fftwf_free(dd);
	fftwf_free(mm);

    	exit(0);
}
