/* Morphological component analysis using 2-D Seislet transform 
Note:  We plan to use analysis based FISTA algorithm. For the time 
being, we simply use ISTA combined with seislet transform. 
Theoreticaly speaking, seislet frame should be better. 	
*/
/*
  Copyright (C) 2004 University of Texas at Austin
   
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

#include "pthresh.h"

int main(int argc, char *argv[])
{
    int niter, n1, n2, nthr; 
    int order;
    float *dobs, *drec1, *drec2, *dtmp, *tmp, *mask;
    float **dip1, **dip2;
    sf_file Fin, Fout1, Fout2, Fmask, Fdip1, Fdip2;
    float pscale, p, pclip, thr, eps;
    char *type, *mode;
    bool unit=false, inv=true, verb;

    sf_init(argc,argv);

    Fin = sf_input("in");
    Fmask=sf_input("mask");  
    Fout1 = sf_output("out");/* component 1 */
    Fout2 = sf_output("comp2");/* component 2 */
    Fdip1=sf_input("dip1");
    Fdip2=sf_input("dip2");

    if (!sf_histint(Fin,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(Fin,"n2",&n2)) sf_error("No n2= in input");

    if(!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */
    if(!sf_getint("order",&order)) order=1;
    /* accuracy order for seislet transform*/
    if (NULL == (type=sf_getstring("type"))) type="linear";
    /* [haar,linear,biorthogonal] wavelet type, the default is linear  */
    if (!sf_getfloat("pscale",&pscale)) pscale=25;
    /* percentile of small scale to be preserved (default is 25)*/

    if(!sf_getbool("verb",&verb))    	verb=false;
    /* verbosity */
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
    if (NULL != sf_getstring("mask")){
	mask=sf_floatalloc(n1*n2);
	sf_floatread(mask,n1*n2,Fmask);
    }

    dobs = sf_floatalloc(n1*n2);
    drec1 = sf_floatalloc(n1*n2);
    drec2 = sf_floatalloc(n1*n2);
    dip1=sf_floatalloc2(n1,n2);
    dip2=sf_floatalloc2(n1,n2);
    dtmp = sf_floatalloc(n1*n2);	
    tmp = sf_floatalloc(n1*n2);

    sf_floatread(dobs,n1*n2,Fin);
    sf_floatread(dip1[0],n1*n2,Fdip1);
    sf_floatread(dip2[0],n1*n2,Fdip2);
    memset(drec1, 0, n1*n2*sizeof(float));
    memset(drec2, 0, n1*n2*sizeof(float));

    seislet_init(n1,n2,inv,unit,eps,order,type[0]);  /* unit=false inv=true */

    for(int iter=1; iter<=niter; iter++)  {
    	seislet_set(dip1);// optimizing component 1 with thresholding 
	for(int i2=0;i2<n2;i2++)
	for(int i1=0; i1<n1; i1++) 
	{	
		float m=(mask[i1+i2*n1])?1:0; 
		drec1[i1+n1*i2]+=dobs[i1+n1*i2]-m*(drec1[i1+n1*i2]+drec2[i1+n1*i2]);
	}
	// seislet adjoint: At(drec)
	seislet_lop(true,false,n1*n2,n1*n2,dtmp,drec1);

	// perform thresholding; T{ At(drec) }
	for(int i2=0;i2<n2;i2++)
	for(int i1=0; i1<n1; i1++) 
	{	
		if (i2>0.01*pscale*n2) dtmp[i1+i2*n1]=0;// set large scale to 0
		tmp[i1+n1*i2]=fabsf(dtmp[i1+n1*i2]);
	}
   	nthr = 0.5+n1*n2*(1.-0.01*pclip);  
    	if (nthr < 0) nthr=0;
    	if (nthr >= n1*n2) nthr=n1*n2-1;
	thr=sf_quantile(nthr,n1*n2,tmp);
	thr*=powf(0.01,(iter-1.0)/(niter-1.0));	// exponentially decrease thr
	pthresholding2(dtmp, n1*n2, thr, p, mode);

	// forward seislet: A T{ At(drec) } 
	seislet_lop(false,false,n1*n2,n1*n2,dtmp,drec1);


    	seislet_set(dip2);// optimizing component 2 with thresholding g 
	for(int i2=0;i2<n2;i2++)
	for(int i1=0; i1<n1; i1++) 
	{	
		float m=(mask[i1+i2*n1])?1:0; 
		drec2[i1+n1*i2]+=dobs[i1+n1*i2]-m*(drec1[i1+n1*i2]+drec2[i1+n1*i2]);
	}
	// seislet adjoint: At(drec)
	seislet_lop(true,false,n1*n2,n1*n2,dtmp,drec2);

	// perform thresholding; T{ At(drec) }
	for(int i2=0;i2<n2;i2++)
	for(int i1=0; i1<n1; i1++) 
	{	
		if (i2>0.01*pscale*n2) dtmp[i1+i2*n1]=0;// set large scale to 0
		tmp[i1+n1*i2]=fabsf(dtmp[i1+n1*i2]);
	}
   	nthr = 0.5+n1*n2*(1.-0.01*pclip);  
    	if (nthr < 0) nthr=0;
    	if (nthr >= n1*n2) nthr=n1*n2-1;
	thr=sf_quantile(nthr,n1*n2,tmp);
	thr*=powf(0.01,(iter-1.0)/(niter-1.0));	// exponentially decrease thr
	pthresholding2(dtmp, n1*n2, thr, p, mode);

	// forward seislet: A T{ At(drec) } 
	seislet_lop(false,false,n1*n2,n1*n2,dtmp,drec2);

	if (verb)    sf_warning("iteration %d;",iter+1);
    }

    sf_floatwrite(drec1,n1*n2,Fout1);
    sf_floatwrite(drec2,n1*n2,Fout2);
    exit(0);
}

