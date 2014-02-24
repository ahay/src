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

#ifdef _OPENMP
#include <omp.h>
#endif


#include "pthresh.h"
#include "ft3d.h"


int main(int argc, char* argv[])
{
    bool verb;
    int niter,n1,n2,n3,i1,i2,i3, iter;
    float p, tol, pclip,t0=1.0, t1, beta, thr;
    char *mode;
    float *din, *mask, *dout;	
    sf_complex *dprev,*dcurr,*dtmp;
    sf_file Fin,Fout, Fmask;

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
    din=sf_floatalloc(n1*n2*n3); 
    dout=sf_floatalloc(n1*n2*n3);
    sf_floatread(din,n1*n2*n3,Fin);
    if (NULL != sf_getstring("mask")){
	mask=sf_floatalloc(n2*n3);
	sf_floatread(mask,n2*n3,Fmask);
    }

    dprev=sf_complexalloc(n1*n2*n3);
    dcurr=sf_complexalloc(n1*n2*n3);
    dtmp=sf_complexalloc(n1*n2*n3);
    for(i1=0; i1<n1*n2*n3; i1++) {
	dprev[i1]=sf_cmplx(din[i1],0.0);
	dcurr[i1]=sf_cmplx(din[i1],0.0);
	dtmp[i1]=sf_cmplx(0.0,0);
    }
    ft3d_init(n1, n2, n3);

    /* FPOCS iterations */
    for(iter=1; iter<=niter; iter++)
    {
	t1=0.5*(1.0+sqrtf(1.0+4.0*t0*t0));
	beta=(t0-1.0)/t1;
	t0=t1;

#ifdef _OPENMP
#pragma omp parallel for default(none)		\
	private(i1)				\
	shared(dtmp,dcurr,beta,dprev,n1,n2,n3)
#endif
	for(i1=0;i1<n1*n2*n3;i1++) {
	    dtmp[i1]=dcurr[i1]+beta*(dcurr[i1]-dprev[i1]);
	    dprev[i1]=dcurr[i1];
	}	


	ft3d_lop(true, false, n1*n2*n3, n1*n2*n3, dcurr, dtmp);

	// perform hard thresholding
	for(i3=0; i3<n3; i3++)
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
	    	dout[i1+n1*i2+n1*n2*i3]=cabsf(dcurr[i1+n1*i2+n1*n2*i3]);
	}

   	int nthr = 0.5+n1*n2*n3*(1.-0.01*pclip); 
    	if (nthr < 0) nthr=0;
    	if (nthr >= n1*n2*n3) nthr=n1*n2*n3-1;
	thr=sf_quantile(nthr,n1*n2*n3,dout);
	thr*=powf(0.01,(iter-1.0)/(niter-1.0));
	sf_cpthresh(dcurr, n1*n2*n3,thr, p,mode);

	ft3d_lop(false, false, n1*n2*n3, n1*n2*n3, dcurr, dtmp);
	
	/* d_rec = d_obs+(1-M)*A T{ At(d_rec) } */
#ifdef _OPENMP
#pragma omp parallel for collapse(3) default(none)	\
	private(i1,i2,i3)				\
	shared(mask,din,dcurr,dtmp,n1,n2,n3)
#endif
	for(i3=0; i3<n3; i3++)	
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{ 
		float m=(mask[i2+i3*n2])?1:0;
		dcurr[i1+n1*(i2+n2*i3)]=sf_cmplx(din[i1+n1*(i2+n2*i3)],0)
			+(1.-m)*dtmp[i1+n1*(i2+n2*i3)];
	}

	if (verb)    sf_warning("iteration %d;",iter);
    }
    for(i1=0;i1<n1*n2*n3; i1++) dout[i1]=crealf(dcurr[i1]);
    sf_floatwrite(dout,n1*n2*n3,Fout);

    ft3d_close();
    sf_close();
    exit(0);
}
