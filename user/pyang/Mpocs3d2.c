/* POCS for 3D missing data interpolation for frequency domain data
NB: using advanced features of FFTW3.
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

void dim_circshift(bool forw, sf_complex *wdat, int n1, int n2, int n3)
/*< circle-shift dimensions for the data:
if forw, shift forward: last-->1st; 1st-->2nd; 2nd-->3rd;...
else, shift backward: 1st-->last; 2nd-->1st; 3rd-->2nd;...
>*/
{
	int i1, i2, i3;
	sf_complex *tmp;
	tmp=(sf_complex*)malloc(n1*n2*n3*sizeof(sf_complex));	
	if(forw){
		for(i3=0; i3<n3; i3++)
		for(i2=0; i2<n2; i2++)
		for(i1=0; i1<n1; i1++)
			tmp[i3+n3*i1+n3*n1*i2]=wdat[i1+n1*i2+n1*n2*i3];
	}else{
		for(i3=0; i3<n3; i3++)
		for(i2=0; i2<n2; i2++)
		for(i1=0; i1<n1; i1++)
			tmp[i2+n2*i3+n2*n3*i1]=wdat[i1+n1*i2+n1*n2*i3];
	}
	memcpy(wdat, tmp, n1*n2*n3*sizeof(sf_complex));
	free(tmp);
}


int main(int argc, char* argv[])
{
    	bool verb;
    	int n1, n2, n3, num, nthr, n[3], i1, i2, i3, index, nw, niter, iter;
    	float thr, pclip, m;		
    	float  *thresh, *tdat, *mask;
	fftwf_complex *wdat, *wdat1;
	fftwf_plan fft1, ifft1, fft23, ifft23;
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
	n[0]=n1;
	n[1]=n2;
	n[2]=n3;	
	nw=n1/2+1;
	num=nw*n2*n3;//total number of elements in frequency domain

    	/* allocate data and mask arrays */
	thresh=(float*)malloc(nw*n2*n3*sizeof(float));
    	tdat=(float*)fftwf_malloc(n1*n2*n3*sizeof(float)); // data in time domain
    	wdat=(fftwf_complex*)fftwf_malloc(nw*n2*n3*sizeof(fftwf_complex));// data in frequency domain
    	wdat1=(fftwf_complex*)fftwf_malloc(nw*n2*n3*sizeof(fftwf_complex));// data in frequency domain
    	fft1=fftwf_plan_many_dft_r2c(1, &n1, n2*n3, tdat, &n1, 1, n1, wdat, &n1, 1, nw, FFTW_MEASURE);	
   	ifft1=fftwf_plan_many_dft_c2r(1, &n1, n2*n3, wdat, &n1, 1, nw, tdat, &n1, 1, n1, FFTW_MEASURE);
	fft23=fftwf_plan_many_dft(2, &n[1], nw, wdat1, &n[1], nw, 1, wdat1, &n[1], nw, 1, FFTW_FORWARD,FFTW_MEASURE);
	ifft23=fftwf_plan_many_dft(2, &n[1], nw, wdat1, &n[1], nw, 1, wdat1, &n[1], nw, 1, FFTW_BACKWARD,FFTW_MEASURE);

	/* initialization */
    	sf_floatread(tdat, n1*n2*n3, Fin);
    	if (NULL != sf_getstring("mask")){
		mask=sf_floatalloc(n2*n3);
		sf_floatread(mask,n2*n3,Fmask);
    	}


	// transform the data from time domain to frequency domain: tdat-->wdat 
	fftwf_execute(fft1);
	for(i1=0; i1<num; i1++) wdat[i1]/=sqrtf(n1);
	memset(wdat1,0,num*sizeof(fftwf_complex));

    	for(iter=1; iter<=niter; iter++)
    	{
		fftwf_execute(fft23);
	#ifdef _OPENMP
	#pragma omp parallel for default(none)	\
		private(i1)			\
		shared(wdat1,num,n2,n3)
	#endif
		for(i1=0; i1<num; i1++) wdat1[i1]/=sqrtf(n2*n3);

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

		fftwf_execute(ifft23);
	#ifdef _OPENMP
	#pragma omp parallel for default(none)	\
		private(i1)			\
		shared(wdat1,num,n2,n3)
	#endif
		for(i1=0; i1<num; i1++) wdat1[i1]/=sqrtf(n2*n3);		

		// d_rec = d_obs+(1-M)*A T{ At(d_rec) } 
	#ifdef _OPENMP
	#pragma omp parallel for collapse(3) default(none)	\
		private(i1,i2,i3,index,m)			\
		shared(mask,wdat,wdat1,nw,n2,n3)
	#endif
		for(i3=0; i3<n3; i3++)	
		for(i2=0; i2<n2; i2++)
		for(i1=0; i1<nw; i1++)
		{ 
			m=(mask[i2+i3*n2])?1:0;
			index=i1+nw*i2+nw*n2*i3;
			wdat1[index]=wdat[index]+(1.-m)*wdat1[index];
		}

		if (verb)    sf_warning("iteration %d;",iter);
    	}

	// transform the data from frequency domain to time domain: wdat-->tdat
	memcpy(wdat, wdat1, num*sizeof(fftwf_complex));
	fftwf_execute(ifft1);
	for(i1=0; i1<n1*n2*n3; i1++) tdat[i1]/=sqrtf(n1);
    	sf_floatwrite(tdat, n1*n2*n3, Fout);

	free(thresh);
	fftwf_free(tdat);
	fftwf_free(wdat);
	fftwf_free(wdat1);		


    	exit(0);
}
