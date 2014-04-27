/* POCS for 5D seismic data interpolation
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
#include <complex.h>
#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    	bool verb;
    	int n1, n2, n3, n4, n5, num, nthr, n[5], i1, i2, i3, i4, i5, index, nw, niter, iter;
    	float thr, pclip, m;		
    	float  *thresh, *tdat, *mask;
	fftwf_complex *wdat, *wdat1;
	fftwf_plan fft1, ifft1, fft2345, ifft2345;
    	sf_file Fin, Fout, Fmask;/* mask and I/O files*/ 

    	sf_init(argc,argv);/* Madagascar initialization */
    	omp_init(); 	/* initialize OpenMP support */

    	/* setup I/O files */
    	Fin=sf_input("in");	/* read the data to be interpolated */
    	Fmask=sf_input("mask");  	/* read the 2-D mask for 3-D data */
    	Fout=sf_output("out"); 	/* output the reconstructed data */
 
    	if(!sf_getbool("verb",&verb)) 	verb=false;
    	/* verbosity */
    	if (!sf_getint("niter",&niter)) 	niter=100;
    	/* total number of POCS iterations */
    	if (!sf_getfloat("pclip",&pclip)) 	pclip=99.;
    	/* starting data clip percentile (default is 99)*/

    	/* Read the data size */
    	if (!sf_histint(Fin,"n1",&n1)) sf_error("No n1= in input");
    	if (!sf_histint(Fin,"n2",&n2)) sf_error("No n2= in input");
    	if (!sf_histint(Fin,"n3",&n3)) sf_error("No n3= in input");
    	if (!sf_histint(Fin,"n4",&n4)) sf_error("No n4= in input");
    	if (!sf_histint(Fin,"n5",&n5)) sf_error("No n5= in input");
	n[0]=n1;
	n[1]=n2;
	n[2]=n3;
	n[3]=n4;
	n[4]=n5;	
	nw=n1/2+1;
	num=nw*n2*n3*n4*n5;//total number of elements in frequency domain

    	/* allocate data and mask arrays */
	thresh=(float*)malloc(nw*n2*n3*n4*n5*sizeof(float));
    	tdat=(float*)fftwf_malloc(n1*n2*n3*n4*n5*sizeof(float)); // data in time domain
    	wdat=(fftwf_complex*)fftwf_malloc(nw*n2*n3*n4*n5*sizeof(fftwf_complex));// data in frequency domain
    	wdat1=(fftwf_complex*)fftwf_malloc(nw*n2*n3*n4*n5*sizeof(fftwf_complex));// data in frequency domain
    	fft1=fftwf_plan_many_dft_r2c(1, &n1, n2*n3*n4*n5, tdat, &n1, 1, n1, wdat, &n1, 1, nw, FFTW_MEASURE);	
   	ifft1=fftwf_plan_many_dft_c2r(1, &n1, n2*n3*n4*n5, wdat, &n1, 1, nw, tdat, &n1, 1, n1, FFTW_MEASURE);
	fft2345=fftwf_plan_many_dft(4, &n[1], nw, wdat1, &n[1], nw, 1, wdat1, &n[1], nw, 1, FFTW_FORWARD,FFTW_MEASURE);
	ifft2345=fftwf_plan_many_dft(4, &n[1], nw, wdat1, &n[1], nw, 1, wdat1, &n[1], nw, 1, FFTW_BACKWARD,FFTW_MEASURE);

	/* initialization */
    	sf_floatread(tdat, n1*n2*n3*n4*n5, Fin);
    	if (NULL != sf_getstring("mask")){
		mask=sf_floatalloc(n2*n3*n4*n5);
		sf_floatread(mask,n2*n3*n4*n5,Fmask);
    	}


	// transform the data from time domain to frequency domain: tdat-->wdat 
	fftwf_execute(fft1);
	for(i1=0; i1<num; i1++) wdat[i1]/=sqrtf(n1);
	memset(wdat1,0,num*sizeof(fftwf_complex));

    	for(iter=1; iter<=niter; iter++)
    	{
		fftwf_execute(fft2345);
	#ifdef _OPENMP
	#pragma omp parallel for default(none)	\
		private(i1)			\
		shared(wdat1,num,n2,n3,n4,n5)
	#endif
		for(i1=0; i1<num; i1++) wdat1[i1]/=sqrtf(n2*n3*n4*n5);

		// perform hard thresholding
	#ifdef _OPENMP
	#pragma omp parallel for default(none)	\
		private(i1)			\
		shared(thresh,wdat1,num)
	#endif
		for(i1=0; i1<num; i1++)	thresh[i1]=cabsf(wdat1[i1]);

	   	nthr = 0.5+num*(1.-0.01*pclip); 
	    	if (nthr < 0) nthr=0;
	    	if (nthr >= num) nthr=num-1;
		thr=sf_quantile(nthr,num,thresh);
		thr*=powf(0.01,(iter-1.0)/(niter-1.0));

	#ifdef _OPENMP
	#pragma omp parallel for default(none)	\
		private(i1)			\
		shared(wdat1,num,thr)
	#endif
		for(i1=0; i1<num; i1++) wdat1[i1]*=(cabsf(wdat1[i1])>thr?1.:0.);

		fftwf_execute(ifft2345);
	#ifdef _OPENMP
	#pragma omp parallel for default(none)	\
		private(i1)			\
		shared(wdat1,num,n2,n3,n4,n5)
	#endif
		for(i1=0; i1<num; i1++) wdat1[i1]/=sqrtf(n2*n3*n4*n5);		

		// d_rec = d_obs+(1-M)*A T{ At(d_rec) } 
	#ifdef _OPENMP
	#pragma omp parallel for collapse(3) default(none)	\
		private(i1,i2,i3,i4,i5,index,m)			\
		shared(mask,wdat,wdat1,nw,n2,n3,n4,n5)
	#endif
		for(i5=0; i5<n5; i5++)	
		for(i4=0; i4<n4; i4++)
		for(i3=0; i3<n3; i3++)	
		for(i2=0; i2<n2; i2++)
		for(i1=0; i1<nw; i1++)
		{ 
			m=(mask[i2+i3*n2+i4*n2*n3+i5*n2*n3*n4])?1:0;
			index=i1+nw*i2+nw*n2*i3+nw*n2*n3*i4+nw*n2*n3*n4*i5;
			wdat1[index]=wdat[index]+(1.-m)*wdat1[index];
		}

		if (verb)    sf_warning("iteration %d;",iter);
    	}

	// transform the data from frequency domain to time domain: wdat-->tdat
	memcpy(wdat, wdat1, num*sizeof(fftwf_complex));
	fftwf_execute(ifft1);
	for(i1=0; i1<n1*n2*n3*n4*n5; i1++) tdat[i1]/=sqrtf(n1);
    	sf_floatwrite(tdat, n1*n2*n3*n4*n5, Fout);

	free(thresh);
	fftwf_free(tdat);
	fftwf_free(wdat);
	fftwf_free(wdat1);		

    	exit(0);
}
