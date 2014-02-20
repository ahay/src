/* 3-D Two-step POCS interpolation using a general Lp-norm optimization
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
    int niter; 
    float p, tol,pclip;
    char *mode;
    sf_file Fin=NULL,Fout=NULL, Fmask=NULL;/* mask and I/O files*/ 

    /* define temporary variables */
    int n1,n2,n3;
    float *din=NULL, *mask=NULL, *dout=NULL;

    sf_init(argc,argv);	/* Madagascar initialization */
#ifdef _OPENMP
    omp_init(); 	/* initialize OpenMP support */
#endif

    /* setup I/O files */
    Fin=sf_input("in");	/* read the data to be interpolated */
    Fmask=sf_input("mask");  	/* read the 2-D mask for 3-D data */
    Fout=sf_output("out"); 	/* output the reconstructed data */
 
    if(!sf_getbool("verb",&verb))    	verb=false;
    /* verbosity */
    if (!sf_getint("niter",&niter)) 	niter=100;
    /* total number iterations */
    if (!sf_getfloat("tol",&tol)) 	tol=1.0e-6;
    /* iteration tolerance */
    if (!sf_getfloat("pclip",&pclip)) 	pclip=99.;
    /* starting data clip percentile (default is 99)*/
    if ( !(mode=sf_getstring("mode")) ) mode = "exp";
    /* thresholding mode: 'hard', 'soft','pthresh','exp';
	'hard', hard thresholding;	'soft', soft thresholding; 
	'pthresh', generalized quasi-p; 'exp', exponential shrinkage */
    if (!sf_getfloat("p",&p)) 		p=0.35;
    /* norm=p, where 0<p<=1 */;
    if (strcmp(mode,"soft") == 0) 	p=1;
    else if (strcmp(mode,"hard") == 0) 	p=0;

    /* Read the data size */
    if (!sf_histint(Fin,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(Fin,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(Fin,"n3",&n3)) sf_error("No n3= in input");

    /* allocate data and mask arrays */
    din=sf_floatalloc(n1*n2*n3); sf_floatread(din,n1*n2*n3,Fin);
    dout=sf_floatalloc(n1*n2*n3);
    if (NULL != sf_getstring("mask")){
	mask=sf_floatalloc(n2*n3);
	sf_floatread(mask,n2*n3,Fmask);
    }


    /********************* 3-D POCS interpolation *********************/
    fftwf_plan p1,p2;/* execute plan for FFT and IFFT */
    int i1,i2,i3, iter;
    float t0=1.0,t1,beta,thr;		
    sf_complex *dprev,*dcurr,*dtmp;

    dprev=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*n1*n2*n3);
    dcurr=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*n1*n2*n3);
    dtmp=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*n1*n2*n3);
    p1=fftwf_plan_dft_3d(n1,n2,n3,dtmp,dtmp,FFTW_FORWARD,FFTW_MEASURE);	
    p2=fftwf_plan_dft_3d(n1,n2,n3,dtmp,dtmp,FFTW_BACKWARD,FFTW_MEASURE);


    for(i1=0; i1<n1*n2*n3; i1++) {
	dprev[i1]=din[i1];
	dcurr[i1]=din[i1];
	dtmp[i1]=0.0;
    }

    /* FPOCS iterations */
    for(iter=1; iter<=niter; iter++)
    {
	t1=0.5*(1.0+sqrtf(1.0+4.0*t0*t0));
	beta=(t0-1.0)/t1;

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i1=0;i1<n1*n2*n3;i1++) {
	    dtmp[i1]=dcurr[i1]+beta*(dcurr[i1]-dprev[i1]);
	    dprev[i1]=dcurr[i1];
	}	

	fftwf_execute(p1);/* FFT */

	/* find the thresholds */
/*
	float mmax=cabsf(dtmp[0]);
	for(i1=1; i1<n1*n2*n3; i1++)
	    mmax=(cabsf(dtmp[i1])>mmax)?cabsf(dtmp[i1]):mmax;
	thr=0.99*powf(0.005,(iter-1.0)/(niter-1.0))*mmax;
*/
	for(i1=1; i1<n1*n2; i1++){
	    dout[i1]=cabsf(dtmp[i1]);
	}

   	int nthr = 0.5+n1*n2*(0.01*pclip);  /*round off*/
    	if (nthr < 0) nthr=0;
    	if (nthr >= n1*n2) nthr=n1*n2-1;
	thr=sf_quantile(nthr,n1*n2,dout);
	thr*=powf(0.01,(iter-1.0)/(niter-1.0));


	/* perform p-norm thresholding */
	pthresholding(dtmp, n1*n2*n3,thr, p,mode);

	fftwf_execute(p2);/* unnormalized IFFT */

	/* adjointness needs scaling with factor 1.0/(n1*n2*n3) */	
#ifdef _OPENMP
#pragma omp parallel for	
#endif
	for(i1=0; i1<n1*n2*n3;i1++) dcurr[i1]=dtmp[i1]/(n1*n2*n3);
	
	/* update d_rec: d_rec = d_obs+(1-M)*A T{ At(d_rec) } */
	for(i3=0;i3<n3;i3++)	
	    for(i2=0;i2<n2;i2++){
		if (mask[i2+i3*n2]){			
		    for(i1=0; i1<n1; i1++)
			dcurr[i1+n1*(i2+n2*i3)]=din[i1+n1*(i2+n2*i3)];
		    }
		}

	if (verb)    sf_warning("iteration %d;",iter);
    }

    /* take the real part */
    for(i1=0;i1<n1*n2*n3; i1++) dout[i1]=crealf(dcurr[i1]);
	
    fftwf_destroy_plan(p1);
    fftwf_destroy_plan(p2);
    fftwf_free(dprev);
    fftwf_free(dcurr);
    fftwf_free(dtmp);

    /* write reconstructed seismic data to output */
    sf_floatwrite(dout,n1*n2*n3,Fout);

    sf_close();
    exit(0);
}
