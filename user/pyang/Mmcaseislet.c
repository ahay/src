/* Morphological component analysis using 2-D Seislet transform 
Note:  We plan to use analysis based FISTA algorithm. For the time 
being, we simply use ISTA combined with seislet transform. 
Theoreticaly speaking, seislet frame should be better. 	
*/

/*
  Copyright (C) 2014 Xi'an Jiaotong University, UT Austin (Pengliang Yang)
   
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
#include <rsfpwd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "pthresh.h"

int main(int argc, char *argv[])
{
    bool verb, smth1, smth2;
    int niter, n1, n2, nthr, order, i1, i2, rect[2], ndat[2];
    float pscale, p, pclip, thr, eps, m;
    float *dobs, *drec1, *drec2, *dtmp, *tmp, *mask;
    float **dip1, **dip2;
    char *type, *mode;
    sf_file Fin, Fout1, Fout2, Fmask, Fdip1, Fdip2;

    sf_init(argc,argv);	/* Madagascar initialization */
#ifdef _OPENMP
    omp_init(); 	/* initialize OpenMP support */
#endif

    Fin = sf_input("in");/* original data */
    Fmask=sf_input("mask");  /* mask for missing values */
    Fout1 = sf_output("out");/* component 1 */
    Fout2 = sf_output("comp2");/* component 2 */
    Fdip1=sf_input("dip1");/* dip of component 1 */
    Fdip2=sf_input("dip2");/*dip of component 2 */

    if (!sf_histint(Fin,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(Fin,"n2",&n2)) sf_error("No n2= in input");

    if(!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */
    if(!sf_getint("order",&order)) order=1;
    /* accuracy order for seislet transform*/
    if (NULL == (type=sf_getstring("type"))) type="linear";
    /* [haar,linear,biorthogonal] wavelet type, the default is linear  */
    if (!sf_getfloat("pscale",&pscale)) pscale=100;
    /* percentile of small scale to be preserved (default is 100)*/

    if(!sf_getbool("verb",&verb))    	verb=false;
    /* verbosity or not */
    if(!sf_getbool("smth1",&smth1))    	smth1=false;
    /* component 1 smoothing or not */
    if(!sf_getbool("smth2",&smth2))    	smth2=false;
    /* component 2 smoothing or not */
    if(!sf_getint("rect1",&rect[0])) 	rect[0]=1;
    /* 1st axis smoothing radius */
    if(!sf_getint("rect2",&rect[1])) 	rect[1]=1;
    /* 2nd axis smoothing radius */
    if (!sf_getint("niter",&niter)) 	niter=30;
    /* total number iterations */
    if (!sf_getfloat("pclip",&pclip)) 	pclip=99;
    /* starting data clip percentile (default is 99)*/
    if ( !(mode=sf_getstring("mode")) ) mode = "exp";
    /* thresholding mode: 'hard', 'soft','pthresh','exp';
	'hard', hard thresholding;	'soft', soft thresholding; 
	'pthresh', generalized quasi-p; 'exp', exponential shrinkage */
    if (!sf_getfloat("p",&p)) 		p=0.35;
    /* norm=p, where 0<p<=1 */;
    if (strcmp(mode,"soft") == 0) 	p=1;
    else if (strcmp(mode,"hard") == 0) 	p=0;

    dobs = sf_floatalloc(n1*n2);
    mask= sf_floatalloc(n1*n2);
    drec1 = sf_floatalloc(n1*n2);
    drec2 = sf_floatalloc(n1*n2);
    dip1= sf_floatalloc2(n1,n2);
    dip2= sf_floatalloc2(n1,n2);
    dtmp = sf_floatalloc(n1*n2);	
    tmp = sf_floatalloc(n1*n2);

    sf_floatread(dobs,n1*n2,Fin);
    sf_floatread(dip1[0],n1*n2,Fdip1);
    sf_floatread(dip2[0],n1*n2,Fdip2);
    memset(drec1, 0, n1*n2*sizeof(float));//memset(drec2, 0, sizeof(*drec1));
    memset(drec2, 0, n1*n2*sizeof(float));//memset(drec2, 0, sizeof(*drec2));
    if (NULL != sf_getstring("mask")){
	sf_floatread(mask,n1*n2,Fmask);
    }else{//no mask, just for separation
	for(int i=0; i<n1*n2; i++) mask[i]=1.;
    }	


    seislet_init(n1,n2,true,false,eps,order,type[0]);//unit=false, inv=true
    if (smth1 || smth2) {//combining shaping with triangle smoothing
	ndat[0]=n1; ndat[1]=n2;
	sf_trianglen_init(2,rect,ndat);
    }

    for(int iter=1; iter<=niter; iter++)  {
	// ========== optimizing component 1 with thresholding =======
    	seislet_set(dip1);

#ifdef _OPENMP
#pragma omp parallel for default(none) collapse(2)	\
	private(i1,i2,m)				\
	shared(n1,n2,drec1,drec2,dobs,mask)
#endif
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++) 
	{	
		m=(mask[i1+i2*n1])?1:0; 
		drec1[i1+n1*i2]=drec1[i1+n1*i2]+dobs[i1+n1*i2]
			-m*(drec1[i1+n1*i2]+drec2[i1+n1*i2]);
	}
	// seislet adjoint: At(drec)
	seislet_lop(true,false,n1*n2,n1*n2,dtmp,drec1);

	// perform thresholding; T{ At(drec) }
#ifdef _OPENMP
#pragma omp parallel for default(none) collapse(2)	\
	private(i1,i2)					\
	shared(n1,n2,dtmp,tmp,pscale)
#endif
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++) 
	{	// set large scale to 0
		if (i2>0.01*pscale*n2) dtmp[i1+i2*n1]=0;
		tmp[i1+n1*i2]=fabsf(dtmp[i1+n1*i2]);
	}
   	nthr = 0.5+n1*n2*(1.-0.01*pclip);  
    	if (nthr < 0) nthr=0;
    	if (nthr >= n1*n2) nthr=n1*n2-1;
	thr=sf_quantile(nthr,n1*n2,tmp);
	thr*=powf(0.01,(iter-1.0)/(niter-1.0));	//exponentially decrease thr
	sf_pthresh(dtmp, n1*n2, thr, p, mode);
	if(smth1){// do smoothing for component 1
		sf_trianglen_lop(true,true,n1*n2,n1*n2,tmp,drec1);
		sf_trianglen_lop(false,false,n1*n2,n1*n2,tmp,drec1);	
	}

	// forward seislet: A T{ At(drec) } 
	seislet_lop(false,false,n1*n2,n1*n2,dtmp,drec1);

	// ============ optimizing component 2 with thresholding =====
    	seislet_set(dip2);
#ifdef _OPENMP
#pragma omp parallel for default(none) collapse(2)	\
	private(i1,i2,m)				\
	shared(n1,n2,drec1,drec2,dobs,mask)
#endif
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++) 
	{	
		m=(mask[i1+i2*n1])?1:0; 
		drec2[i1+n1*i2]=drec2[i1+n1*i2]+dobs[i1+n1*i2]
			-m*(drec1[i1+n1*i2]+drec2[i1+n1*i2]);
	}
	// seislet adjoint: At(drec)
	seislet_lop(true,false,n1*n2,n1*n2,dtmp,drec2);

	// perform thresholding; T{ At(drec) }
#ifdef _OPENMP
#pragma omp parallel for default(none) collapse(2)	\
	private(i1,i2)					\
	shared(n1,n2,dtmp,tmp,pscale)
#endif
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++) 
	{	// set large scale to 0
		if (i2>0.01*pscale*n2) dtmp[i1+i2*n1]=0;
		tmp[i1+n1*i2]=fabsf(dtmp[i1+n1*i2]);
	}
   	nthr = 0.5+n1*n2*(1.-0.01*pclip);  
    	if (nthr < 0) nthr=0;
    	if (nthr >= n1*n2) nthr=n1*n2-1;
	thr=sf_quantile(nthr,n1*n2,tmp);
	thr*=powf(0.01,(iter-1.0)/(niter-1.0));	//exponentially decrease thr
	sf_pthresh(dtmp, n1*n2, thr, p, mode);
	if(smth2){// do smoothing for component 2
		sf_trianglen_lop(true,true,n1*n2,n1*n2,tmp,drec2);
		sf_trianglen_lop(false,false,n1*n2,n1*n2,tmp,drec2);	
	}

	// forward seislet: A T{ At(drec) } 
	seislet_lop(false,false,n1*n2,n1*n2,dtmp,drec2);

	if (verb)    sf_warning("iteration %d;",iter);
    }

    sf_floatwrite(drec1,n1*n2,Fout1);
    sf_floatwrite(drec2,n1*n2,Fout2);

    sf_trianglen_close();
    exit(0);
}

