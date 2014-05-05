/* Morphological component analysis using 2-D Seislet transform 
Note:  We plan to use analysis based iterative shrinkage-thresholding (IST)
 algorithm. Here, nc components with nc seislet transform builds a seislet 
 frame to do the simultineous multicomponent separation and interpolation.	
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

References: 
	1) Elad, Michael, et al. "Simultaneous cartoon and texture image 
	inpainting using morphological component analysis (MCA)." Applied 
	and Computational Harmonic Analysis 19.3 (2005): 340-358.
	2) Starck, Jean-Luc, Michael Elad, and David Donoho. "Redundant 
	multiscale transforms and their application for morphological 
	component separation." Advances in Imaging and Electron Physics
	132.82 (2004): 287-348.

To know why MCA algorithm work like this, it will be much easier if you see 
' Yang, Pengliang, Jinghuai Gao, and Wenchao Chen. "L1/2-constrained 
morphological component analysis." 2013 IEEE China Summit & International 
Conference on Signal and Information Processing (ChinaSIP), 2013. '. The
only difference lies in the thresholding function and the transform used. 
*/

#include <rsf.h>
#include <rsfpwd.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "pthresh.h"

int main(int argc, char *argv[])
{
    bool verb, decr;
    int niter, n1, n2, nc, nthr, order, iter, i1, i2, ic;
    float pscale, p, pclip, thr, eps, m;
    float *dobs, *drec, *mask, *dips, *coeffs, *res, *tmp;
    float **dip;
    char *type, *mode;
    sf_file Fin, Fout, Fmask, Fdips;

    sf_init(argc,argv);	/* Madagascar initialization */

    Fin = sf_input("in");/* original data */
    Fout = sf_output("out");/* nc components */
    Fdips=sf_input("dips");/* dips of nc component */

    if (!sf_histint(Fin,"n1",&n1)) 	sf_error("No n1= in input");
    if (!sf_histint(Fin,"n2",&n2)) 	sf_error("No n2= in input");
    if (!sf_histint(Fdips,"n3",&nc)) 	nc=1;/* nc dips */
    if(!sf_getfloat("eps",&eps)) 	eps=0.01;
    /* regularization */
    if(!sf_getint("order",&order)) 	order=1;
    /* accuracy order for seislet transform*/
    if (NULL == (type=sf_getstring("type"))) type="linear";
    /* [haar,linear,biorthogonal] wavelet type, the default is linear  */
    if (!sf_getfloat("pscale",&pscale)) pscale=100;
    /* percentile of small scale to be preserved (default is 100)*/
    if(!sf_getbool("verb",&verb))    	verb=false;
    /* verbosity or not */
    if(!sf_getbool("decr",&decr))    	decr=true;
    /* decrease threshold in iterations or not */
    if (!sf_getint("niter",&niter)) 	niter=10;
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

    sf_putint(Fout,"n1",n1);
    sf_putint(Fout,"n2",n2);
    sf_putint(Fout,"n3",nc);

    dobs=(float*)malloc(n1*n2*sizeof(float));	// observations
    res=(float*)malloc(n1*n2*sizeof(float));	// residual 
    drec=(float*)malloc(nc*n1*n2*sizeof(float));// reconstructed nc components
    dips=(float*)malloc(nc*n1*n2*sizeof(float));// all estimated dips
    dip=sf_floatalloc2(n1, n2);			// current dip
    coeffs=(float*)malloc(n1*n2*sizeof(float));	// seislet coefficients
    tmp=(float*)malloc(n1*n2*sizeof(float));	// store the absolute
    mask=(float*)malloc(n1*n2*sizeof(float));	// mask

    sf_floatread(dobs, n1*n2, Fin);
    sf_floatread(dips, n1*n2*nc, Fdips);
    memset(drec, 0, n1*n2*nc*sizeof(float));
    memset(coeffs, 0, n1*n2*sizeof(float));
    memset(tmp, 0, n1*n2*sizeof(float));
    memset(res, 0, n1*n2*sizeof(float));
    memset(dip[0], 0, n1*n2*sizeof(float));

    if (NULL != sf_getstring("mask")){
    	Fmask=sf_input("mask");  /* mask for missing values */
	sf_floatread(mask, n1*n2, Fmask);
    }else{//no mask, just for separation
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++) 
	mask[i1+i2*n1]=1;
    }	


    seislet_init(n1, n2, true, false, eps, order, type[0]);//unit=false, inv=true
    seislet_set(dip);
    for(iter=0; iter<niter; iter++)
    {
	memset(res, 0, n1*n2*sizeof(float));
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++) 
	{	
		//sum of nc components-->res
		for(ic=0; ic<nc; ic++) res[i1+n1*i2]+=drec[i1+n1*i2+n1*n2*ic];
		m=(mask[i1+i2*n1])?1.:0; 
		res[i1+n1*i2]=dobs[i1+n1*i2]-m*res[i1+n1*i2];
	}

	for(ic=0; ic<nc; ic++)
	{
		// sparsifying ic-th component with shrinkage/thresholding
		memcpy(dip[0], &dips[ic*n1*n2], n1*n2*sizeof(float));

		for(i2=0; i2<n2; i2++)
		for(i1=0; i1<n1; i1++) 
		{	
			drec[i1+n1*i2+n1*n2*ic]+=res[i1+n1*i2];
		}
		// seislet adjoint: At(drec^{ic})
		seislet_lop(true, false, n1*n2, n1*n2, coeffs, &drec[ic*n1*n2]);

		// perform thresholding; T{ At(drec) }
		for(i2=0; i2<n2; i2++)
		for(i1=0; i1<n1; i1++) 
		{
			if (i2>0.01*pscale*n2) coeffs[i1+i2*n1]=0;// set large scale to 0
			tmp[i1+n1*i2]=fabsf(coeffs[i1+n1*i2]);
		}
		nthr = 0.5+n1*n2*(1.-0.01*pclip);  
		if (nthr < 0) nthr=0;
		if (nthr >= n1*n2) nthr=n1*n2-1;
		thr=sf_quantile(nthr, n1*n2, tmp);
		if(decr) thr*=(float)(niter-iter)/niter;
		sf_pthresh(coeffs, n1*n2, thr, p, mode);

		// forward seislet: A T{ At(drec^{ic}) } 
		seislet_lop(false, false, n1*n2, n1*n2, coeffs, &drec[ic*n1*n2]);
	}

	if (verb) sf_warning("iteration %d;",iter+1);	
    }

    sf_floatwrite(drec, n1*n2*nc, Fout);

    free(dobs);
    free(res);
    free(drec);
    free(dips);
    free(*dip); free(dip);
    free(coeffs);
    free(tmp);
    free(mask);

    exit(0);
}

