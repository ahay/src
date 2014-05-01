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
#include "fftn.h"


int main(int argc, char* argv[])
{
    bool verb;
    int n1,n2,n3,i1,i2,i3,index, num, nthr, niter, iter;
    int n[3];
    float p, tol, pclip, t0, t1, beta, thr,m;
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

    num=n1*n2*n3;
    n[0]=n1; n[1]=n2; n[2]=n3;

    /* allocate data and mask arrays */
    din=sf_floatalloc(num); 
    dout=sf_floatalloc(num);
    sf_floatread(din,num,Fin);
    if (NULL != sf_getstring("mask")){
	mask=sf_floatalloc(n2*n3);
	sf_floatread(mask,n2*n3,Fmask);
    }

    dprev=sf_complexalloc(num);
    dcurr=sf_complexalloc(num);
    dtmp=sf_complexalloc(num);
    for(i1=0; i1<num; i1++) {
	dprev[i1]=sf_cmplx(din[i1],0.0);
	dcurr[i1]=sf_cmplx(din[i1],0.0);
	dtmp[i1]=sf_cmplx(0.0,0);
    }
    fftn_init(3, n);

    /* FPOCS iterations */
    t0=1.;
    for(iter=0; iter<niter; iter++)
    {
	t1=0.5*(1.0+sqrtf(1.0+4.0*t0*t0));
	beta=(t0-1.0)/t1;
	t0=t1;

#ifdef _OPENMP
#pragma omp parallel for default(none)		\
	private(i1)				\
	shared(dtmp,dcurr,beta,dprev,num)
#endif
	for(i1=0;i1<num;i1++) {
	    dtmp[i1]=dcurr[i1]+beta*(dcurr[i1]-dprev[i1]);
	    dprev[i1]=dcurr[i1];
	}	

	fftn_lop(true, false, num, num, dcurr, dtmp);

	// perform hard thresholding
#ifdef _OPENMP
#pragma omp parallel for default(none)	\
	private(i1)			\
	shared(dout,dcurr,num)
#endif
	for(i1=0; i1<num; i1++)	dout[i1]=cabsf(dcurr[i1]);

   	nthr = 0.5+num*(1.-0.01*pclip); 
    	if (nthr < 0) nthr=0;
    	if (nthr >= num) nthr=num-1;
	thr=sf_quantile(nthr,num,dout);
	//thr*=powf(0.01, iter/(niter-1));
	sf_cpthresh(dcurr, num,thr, p,mode);

	fftn_lop(false, false, num, num, dcurr, dtmp);
	
	/* d_rec = d_obs+(1-M)*A T{ At(d_rec) } */
#ifdef _OPENMP
#pragma omp parallel for collapse(3) default(none)	\
	private(i1,i2,i3,index,m)				\
	shared(mask,din,dcurr,dtmp,n1,n2,n3)
#endif
	for(i3=0; i3<n3; i3++)	
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{ 
		m=(mask[i2+i3*n2])?1:0;
		index=i1+n1*i2+n1*n2*i3;
		dcurr[index]=din[index]	+(1.-m)*dtmp[index];
	}

	if (verb)    sf_warning("iteration %d;",iter+1);
    }

    for(i1=0; i1<num; i1++) dout[i1]=crealf(dcurr[i1]);
    sf_floatwrite(dout,num,Fout);

    fftn_close();
    sf_close();
    exit(0);
}
