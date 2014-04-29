/* Linear/parabolic radon transform frequency domain implementation and its adjoint
Also referred to as high-resolution radon transform
Note: I borrowed a lot from /system/seismic/radon+Mradon.c. The distinction:
	I am using FFTW because I am inexperienced in invoking kiss_fft. 
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

  References: 
	1) Kostov C., 1990. "Toeplitz structure in Slant-Stack inversion":
	 SEG Extended Abstracts, 1647-1650.
	2) Sacchi, Mauricio D., and Milton Porsani. "Fast high resolution 
	parabolic Radon transform." Society of Exploration Geophysicists 
	69th Annual International Meeting, SPRO P. Vol. 1. No. 1. 1999.
	3) Vogel, Curtis R. Computational methods for inverse problems. 
	Vol. 23. Siam, 2002.	[Chapter 5.2.]
*/

#include <rsf.h>
#include <time.h>
#include <complex.h>
#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "myradon2.h"

void matrix_transpose(sf_complex *matrix, int nx, int nz)
{
	int ix, iz;
	sf_complex *tmp=(sf_complex*)malloc(nx*nz*sizeof(sf_complex));
	if (tmp==NULL) {printf("out of memory!\n"); exit(1);}
	for(iz=0; iz<nz; iz++)
	for(ix=0; ix<nx; ix++)
		tmp[iz+nz*ix]=matrix[ix+nx*iz];

	memcpy(matrix, tmp, nx*nz*sizeof(sf_complex));
	free(tmp);
}

int main(int argc, char* argv[])
{
	bool adj, inv, par;
    	char *invmode;
	int iw, ip, ix, np, nt, nfft, nw, nx, niter;
	float dp, p0, dt, t0, dx, ox, x0, w, eps;
	float *p, *xx, **dd, **mm, *tmpr;
	sf_complex *cdd, *cmm;
	fftwf_complex *tmpc;
	fftwf_plan fft1, ifft1;
	sf_file in, out;

    	sf_init(argc,argv);
	in = sf_input("in");	/* veloctiy model */
	out =sf_output("out");	/* shot records */

    	if (!sf_getbool("adj",&adj)) adj=true;
	/* if y, perform adjoint operation */
    	if (!sf_getbool("inv",&inv)) inv=adj; 
	/* if y, perform inverse operation */

    	/* read input file parameters */
    	if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
	/* number of samples in time axis */
    	if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
	/* interval of time axis */
    	if (!sf_histfloat(in,"o1",&t0)) t0=0.;
	/* origin of time axis */
	nfft=2*kiss_fft_next_fast_size(nt);
	nw=nfft/2+1;

    	if (adj||inv) { // m(tau,p)=sum_{i=0}^{nx} d(t=tau+p*x_i,x_i)
		if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
		/* number of offset if the input in the data domain */

		/* specify slope axis */
		if (!sf_getint  ("np",&np)) sf_error("Need np=");
		/* number of p values (if adj=y) */
		if (!sf_getfloat("dp",&dp)) sf_error("Need dp=");
		/* p sampling (if adj=y) */
		if (!sf_getfloat("p0",&p0)) sf_error("Need p0=");
		/* p origin (if adj=y) */
		if(inv){
			
    			if ( !(invmode=sf_getstring("invmode")) ) invmode="toeplitz";
			/* inverse method: 'ls' if least-squares; 'toeplitz' if use FFT */			
			if ((invmode[0]=='l') && !sf_getint("niter",&niter)) niter=100;
			/* number of CGLS iterations */
			if (!sf_getfloat("eps",&eps)) eps=0.01;
			/* regularization parameter */
		}

		sf_putint(  out,"n2",np);
		sf_putfloat(out,"d2",dp);
		sf_putfloat(out,"o2",p0);
    	} else { // d(t,h)=sum_{i=0}^{np} m(tau=t-p_i*h,p_i)
		if (!sf_histint  (in,"n2",&np)) sf_error("No n2= in input");
		/* number of ray parameter if input in radon domain */
		if (!sf_histfloat(in,"d2",&dp)) sf_error("No d2= in input");
		/* p sampling interval if input in radon domain */
		if (!sf_histfloat(in,"o2",&p0)) sf_error("No o2= in input");
		/* p origin if input in radon domain */
		if (!sf_getint("nx",&nx)) sf_error ("Need nx=");
		/* number of offsets (if adj=n) */
	
		sf_putint(out,"n2",nx);
    	}

	p=sf_floatalloc(np);
	xx=sf_floatalloc(nx);
	dd=sf_floatalloc2(nt, nx);
	mm=sf_floatalloc2(nt, np);
	cdd=(sf_complex*)malloc(nw*nx*sizeof(sf_complex));
	cmm=(sf_complex*)malloc(nw*np*sizeof(sf_complex));
	tmpr=(float*)fftwf_malloc(nfft*sizeof(float));
    	tmpc=(fftwf_complex*)fftwf_malloc(nw*sizeof(fftwf_complex));
    	fft1=fftwf_plan_dft_r2c_1d(nfft,tmpr,tmpc,FFTW_MEASURE);	
   	ifft1=fftwf_plan_dft_c2r_1d(nfft,tmpc,tmpr,FFTW_MEASURE);

	for(ip=0; ip<np; ip++) p[ip]=p0+ip*dp;	
	if (adj) {// m(tau,p)=sum_{i=0}^{nx} d(t=tau+p*x_i,x_i)
		sf_floatread(dd[0], nt*nx, in);

	    	if (!sf_histfloat(in,"o2",&ox)) sf_error("No o2= in input");
		/* data origin in x */
	    	if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
		/* sampling interval in x */
	} else {// d(t,h)=sum_{i=0}^{np} m(tau=t-p_i*h,p_i)
		sf_floatread(mm[0], nt*np, in);
	    	if (!sf_getfloat("ox",&ox)) sf_error("Need ox=");
		/* x origin */
	    	if (!sf_getfloat("dx",&dx)) sf_error("Need dx=");
		/* sampling interval in x */

	    	sf_putfloat(out,"o2",ox);
	    	sf_putfloat(out,"d2",dx);
	}
	for(ix=0; ix<nx; ix++) xx[ix] = ox + ix*dx;

    	if (!sf_getbool("parab",&par)) par=false;
	/* if y, parabolic Radon transform */
    	if (!sf_getfloat("x0",&x0)) x0=1.;   
	/* reference offset */

	for (ix=0; ix < nx; ix++)/* normalize offsets */
	{
		if (par) xx[ix] *= xx[ix]/(x0*x0);
		else if (x0!=1.) xx[ix] /= x0;
	}

	if(adj||inv){// m(tau,p)=sum_{i=0}^{nx} d(t=tau+p*x_i,x_i)
		for(ix=0; ix<nx; ix++) // loop over offsets
		{
			memset(tmpr, 0, nfft*sizeof(float));
			memcpy(tmpr, dd[ix], nt*sizeof(float));
		 	fftwf_execute(fft1);// FFT: dd-->cdd
			memcpy(&cdd[ix*nw], tmpc, nw*sizeof(sf_complex));
		}
		matrix_transpose(cdd, nw, nx);
	}else{	// d(t,h)=sum_{i=0}^{np} m(tau=t-p_i*h,p_i)
		for(ip=0; ip<np; ip++) // loop over slopes
		{
			memset(tmpr, 0, nfft*sizeof(float));
			memcpy(tmpr, mm[ip], nt*sizeof(float));
		 	fftwf_execute(fft1);// FFT: mm-->cmm
			memcpy(&cmm[ip*nw], tmpc, nw*sizeof(float));			
		}
		matrix_transpose(cmm, nw, np);
	}


	myradon2_init(np, nx, dp, p, xx);
	for(iw=0; iw<nw; iw++) 
	{
		w=2.*SF_PI*iw/(nfft*dt);
		myradon2_set(w);
		myradon2_lop(adj, false, np, nx, &cmm[iw*np], &cdd[iw*nx]);
		if(inv){
			if (invmode[0]=='t') 
				myradon2_inv(&cmm[iw*np], &cmm[iw*np], eps);
			else
				sf_csolver_reg(myradon2_lop, sf_ccgstep, sf_ccopy_lop, np, 
				np, nx, &cmm[iw*np], &cdd[iw*nx], niter,eps,"end");
		}
	}


	if(adj||inv){// m(tau,p)=sum_{i=0}^{nx} d(t=tau+p*x_i,x_i)
		matrix_transpose(cmm, np, nw);
		for(ip=0; ip<np; ip++) // loop over slopes
		{			
			memcpy(tmpc, &cmm[ip*nw], nw*sizeof(sf_complex));
		 	fftwf_execute(ifft1); // IFFT: cmm-->mm
			//memcpy(mm[ip], tmpr, nt*sizeof(float));
			for(iw=0; iw<nt; iw++) mm[ip][iw]=tmpr[iw]/nfft;
		}

		sf_floatwrite(mm[0], nt*np, out);
	}else{// d(t,h)=sum_{i=0}^{np} m(tau=t-p_i*h,p_i)
		matrix_transpose(cdd, nx, nw);
		for(ix=0; ix<nx; ix++) // loop over offsets
		{
			memcpy(tmpc, &cdd[ix*nw], nw*sizeof(sf_complex));
		 	fftwf_execute(ifft1);// IFFT: cmm-->mm
			//memcpy(dd[ix], tmpr, nt*sizeof(float));
			for(iw=0; iw<nt; iw++) dd[ix][iw]=tmpr[iw]/nfft;
		}

		sf_floatwrite(dd[0], nt*nx, out);
	}


	free(p);
	free(xx);
	free(*dd); free(dd);
	free(*mm); free(mm);
	free(cdd); 
	free(cmm); 
	fftwf_free(tmpr);
	fftwf_free(tmpc);
	fftwf_destroy_plan(fft1);
    	fftwf_destroy_plan(ifft1);

    	exit(0);
}

