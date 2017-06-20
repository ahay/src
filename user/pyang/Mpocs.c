/* n-D POCS interpolation using a hard thresholding
Note: Acquistion geometry specified by mask operator.
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
#include <math.h>
#include <complex.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef SF_HAS_FFTW
#include <fftw3.h>

int main(int argc, char* argv[])
{
    	bool verb;
    	int i, i1, i2, index, n1, n2, num, dim, n[SF_MAX_DIM], nw, iter, niter, nthr;
    	float thr, pclip;
    	float *dobs_t, *thresh, *mask;
    	char key[7];
    	fftwf_complex *mm, *dd, *dobs;
    	fftwf_plan fft1, ifft1, fft2, ifft2;/* execute plan for FFT and IFFT */
    	sf_file in, out, Fmask;/* mask and I/O files*/ 


    	sf_init(argc,argv);	/* Madagascar initialization */
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
	num=nw*n2;/* total number of elements in frequency domain */
 
    	/* allocate data and mask arrays */
	thresh=(float*)malloc(nw*n2*sizeof(float));
    	dobs_t=(float*)fftwf_malloc(n1*n2*sizeof(float)); /* data in time domain*/
    	dobs=(fftwf_complex*)fftwf_malloc(nw*n2*sizeof(fftwf_complex));
    	dd=(fftwf_complex*)fftwf_malloc(nw*n2*sizeof(fftwf_complex));
    	mm=(fftwf_complex*)fftwf_malloc(nw*n2*sizeof(fftwf_complex));

    	fft1=fftwf_plan_many_dft_r2c(1, &n1, n2, dobs_t, &n1, 1, n1, dobs, &n1, 1, nw, FFTW_MEASURE);	
   	ifft1=fftwf_plan_many_dft_c2r(1, &n1, n2, dobs, &n1, 1, nw, dobs_t, &n1, 1, n1, FFTW_MEASURE);
	fft2=fftwf_plan_many_dft(dim-1, &n[1], nw, mm, &n[1], nw, 1, mm, &n[1], nw, 1, FFTW_FORWARD,FFTW_MEASURE);
	ifft2=fftwf_plan_many_dft(dim-1, &n[1], nw, mm, &n[1], nw, 1, mm, &n[1], nw, 1, FFTW_BACKWARD,FFTW_MEASURE);

	/* initialization */
    	sf_floatread(dobs_t, n1*n2, in);
    	if (NULL != sf_getstring("mask")){
		mask=sf_floatalloc(n2);
		sf_floatread(mask,n2,Fmask);
    	} else {
	    mask = NULL;
	    sf_error("Need mask=");
	}

	/*transform the data from time domain to frequency domain: tdat-->wdat*/
	fftwf_execute(fft1);
	for(i=0; i<num; i++) dobs[i]/=sqrtf(n1);
	memset(dd,0,num*sizeof(fftwf_complex));

	/* Projection onto convex sets (POCS) Algorithm:dd^{k+1}=dobs+(1-M)AT[A^* dd^k]	*/
    	for(iter=0; iter<niter; iter++)
    	{
		/* mm<--A^t dd */
		memcpy(mm, dd, num*sizeof(fftwf_complex));
		fftwf_execute(fft2);
		for(i=0; i<num; i++) mm[i]/=sqrtf(n2);

		/* perform hard thresholding: mm<--T{mm} */
		for(i=0; i<num; i++)	thresh[i]=cabsf(mm[i]);
	   	nthr = 0.5+num*(1.-0.01*pclip); 
	    	if (nthr < 0) nthr=0;
	    	if (nthr >= num) nthr=num-1;
		thr=sf_quantile(nthr,num,thresh);
		for(i=0; i<num; i++) mm[i]*=(cabsf(mm[i])>thr?1.:0.);

		/* mm<--A mm*/
		fftwf_execute(ifft2);
		for(i=0; i<num; i++) mm[i]/=sqrtf(n2);		

		/* dd^{k+1}=dobs+(1-M)AT[A^* dd^k] */
		for(i2=0; i2<n2; i2++)
		for(i1=0; i1<nw; i1++)
		{ 
			index=i1+nw*i2;
			dd[index]=dobs[index]+(1.-mask[i2])*mm[index];
		}

		if (verb)    sf_warning("iteration %d;",iter+1);
    	}

	/*transform the data from frequency domain to time domain: wdat-->tdat*/
	memcpy(dobs, dd, num*sizeof(fftwf_complex));
	fftwf_execute(ifft1);
	for(i=0; i<n1*n2; i++) dobs_t[i]/=sqrtf(n1);
    	sf_floatwrite(dobs_t, n1*n2, out);/* output reconstructed seismograms */

	free(thresh);
	fftwf_free(dobs_t);	
	fftwf_free(dobs);
	fftwf_free(mm);
	fftwf_free(dd);

    	exit(0);
}
#endif
