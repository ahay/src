/* Morphological component analysis using 2-D Radon transform	
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
#include "radonlop.h"

int main(int argc, char *argv[])
{
	bool verb, decr,par=true;
	int ip, ix, np, nt, nx, niter, nthr, iter, i1, i2, ic;
	float dp, p0, dt, t0, dx, ox, x0, eps, normp,pclip, thr, m;
	float *p, *xx, *dobs, *drec, *mask, *coeffs, *res, *tmp;
    	char *mode;
	sf_file in, out, offset, Fmask;

    	sf_init(argc,argv);
	in = sf_input("in");	/* whole data */
	out =sf_output("out");	/* separated linear and parabolic components */

    	/* read input file parameters */
    	if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
	/* number of samples in time axis */
    	if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
	/* interval of time axis */
    	if (!sf_histfloat(in,"o1",&t0)) t0=0.;
	/* origin of time axis */

	if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
	/* number of offset if the input in the data domain */
    	if (!sf_histfloat(in,"o2",&ox)) sf_error("No o2= in input");
	/* data origin in x */
    	if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
	/* sampling interval in x */
    	if (!sf_getfloat("x0",&x0)) x0=1.;   
	/* reference offset */

	/* specify slope axis */
	if (!sf_getint  ("np",&np)) sf_error("Need np=");
	/* number of p values (if adj=y) */
	if (!sf_getfloat("dp",&dp)) sf_error("Need dp=");
	/* p sampling (if adj=y) */
	if (!sf_getfloat("p0",&p0)) sf_error("Need p0=");
	/* p origin (if adj=y) */
	if (!sf_getfloat("eps",&eps)) eps=0.01;
	/* regularization parameter */
    	if (!sf_getbool("parab",&par)) par=false;
	/* if y, parabolic Radon transform */

    	if(!sf_getbool("verb",&verb))    	verb=false;
    	/* verbosity or not */
    	if(!sf_getbool("decr",&decr))    	decr=true;
    	/* decrease threshold in iterations or not */
    	if (!sf_getint("niter",&niter)) 	niter=20;
    	/* total number iterations */
    	if (!sf_getfloat("pclip",&pclip)) 	pclip=99;
    	/* starting data clip percentile (default is 99)*/
    	if ( !(mode=sf_getstring("mode")) ) 	mode = "exp";
    	/* thresholding mode: 'hard', 'soft','pthresh','exp';
	'hard', hard thresholding;	'soft', soft thresholding; 
	'pthresh', generalized quasi-p; 'exp', exponential shrinkage */
    	if (!sf_getfloat("p",&normp)) 		normp=0.35;
    	/* norm=p, where 0<p<=1 */;
   	if (strcmp(mode,"soft") == 0) 		normp=1;
    	else if (strcmp(mode,"hard") == 0) 	normp=0;

    	if (NULL != sf_getstring("offset")) {
		offset = sf_input("offset");
		if (nx!=sf_filesize(offset)) sf_error("Wrong dimensions in offset");
    	} else {
		offset = NULL;
    	}

    	sf_putint(out,"n1",nt);
	sf_putfloat(out, "o1", t0);
	sf_putfloat(out, "d1", dt);
    	sf_putint(out,"n2",nx);
	sf_putfloat(out, "o2", ox);
	sf_putfloat(out, "d2", dx);
    	sf_putint(out,"n3",2);

	p=(float*)malloc(np*sizeof(float));		// ray parameter
	xx=(float*)malloc(nx*sizeof(float));		// offset
    	dobs=(float*)malloc(nt*nx*sizeof(float));	// observations
    	res=(float*)malloc(nt*nx*sizeof(float));	// residual 
    	drec=(float*)malloc(2*nt*nx*sizeof(float));	// reconstructed nc components
    	coeffs=(float*)malloc(nt*np*sizeof(float));	// seislet coefficients
    	tmp=(float*)malloc(nt*np*sizeof(float));	// absolute valuse of coeffs
    	mask=(float*)malloc(nt*nx*sizeof(float));	// mask

	for(ip=0; ip<np; ip++) p[ip]=p0+ip*dp;	
    	if (NULL != offset) {
		sf_floatread(xx,nx,offset);
		sf_fileclose(offset);
    	} else {
		for(ix=0; ix<nx; ix++) xx[ix]=ox+ix*dx;
	}
    	sf_floatread(dobs, nt*nx, in);
    	memset(drec, 0, 2*nt*nx*sizeof(float));
    	memset(coeffs, 0, nt*np*sizeof(float));
    	memset(tmp, 0, nt*np*sizeof(float));
    	memset(res, 0, nt*nx*sizeof(float));
    	if (NULL != sf_getstring("mask")){
	    	Fmask=sf_input("mask");  /* mask for missing values */
		sf_floatread(mask, nt*nx, Fmask);
    	}else{//no mask, just for separation
		for(i2=0; i2<nx; i2++)
		for(i1=0; i1<nt; i1++) 
			mask[i1+i2*nt]=1;
    	}	
	sf_radon2_init(true, p, xx, nt, np, nx, x0, dt, dp, eps);

    	for(iter=0; iter<niter; iter++)
    	{
		memset(res, 0, nt*nx*sizeof(float));
		for(i2=0; i2<nx; i2++)
		for(i1=0; i1<nt; i1++) 
		{	
			//sum of 2 components-->res
			for(ic=0; ic<2; ic++) res[i1+nt*i2]+=drec[i1+nt*i2+nt*nx*ic];
			m=(mask[i1+i2*nt])?1.:0; 
			res[i1+nt*i2]=dobs[i1+nt*i2]-m*res[i1+nt*i2];
		}

		for(ic=0; ic<2; ic++)
		{
			sf_radon2_set(par);
			// sparsifying ic-th component with shrinkage/thresholding
			for(i2=0; i2<nx; i2++)
			for(i1=0; i1<nt; i1++) 
			{	
				drec[i1+nt*i2+nt*nx*ic]+=res[i1+nt*i2];
			}
			// seislet adjoint: At(drec^{ic})
			sf_radon2_lop(true, false, nt*np, nt*nx, coeffs, &drec[ic*nt*nx]);

			// perform thresholding; T{ At(drec) }
			for(i2=0; i2<np; i2++)
			for(i1=0; i1<nt; i1++) 
				tmp[i1+nt*i2]=fabsf(coeffs[i1+nt*i2]);
			nthr = 0.5+nt*np*(1.-0.01*pclip);  
			if (nthr < 0) nthr=0;
			if (nthr >= nt*np) nthr=nt*np-1;
			thr=sf_quantile(nthr, nt*np, tmp);
			sf_warning("thr=%g",thr);
			if(decr) thr*=(float)(niter-iter)/niter;
			sf_warning("thr=%g",thr);
			sf_pthresh(coeffs, nt*np, thr, normp, mode);

			// forward seislet: A T{ At(drec^{ic}) } 		
			sf_radon2_lop(false, false, nt*np, nt*nx, coeffs, &drec[ic*nt*nx]);

			par=!par;//parabolic-->linear; linear-->parabolic
		}

		if (verb) sf_warning("iteration %d;",iter+1);	
    	}

    	sf_floatwrite(drec, nt*nx*2, out);

	free(p);
	free(xx);
    	free(dobs);
    	free(res);
    	free(drec);
    	free(coeffs);
    	free(tmp);
    	free(mask);

    	exit(0);
}

