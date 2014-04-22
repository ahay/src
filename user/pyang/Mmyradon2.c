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

  Reference: Kostov C., 1990. "Toeplitz structure in Slant-Stack
	inversion": SEG Extended Abstracts, 1647-1650.
*/

#include <rsf.h>
#include <time.h>
#include <complex.h>
#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
	bool adj, inv, verb, par;
	int it, ip, ix, np, nt, nfft, nx;
	float dp, p0, dt, t0, dx, ox, x0;
	float *p, *xx, **dd, **mm;
	fftwf_complex *tmp;
	fftwf_plan fft1, ifft1;
	sf_complex tmpc, **cdd, **cmm;
	sf_file in, out;

    	sf_init(argc,argv);
	in = sf_input("in");	/* veloctiy model */
	out =sf_output("out");	/* shot records */

    	if (!sf_getbool("adj",&adj)) adj=true;
	/* if y, perform adjoint operation */
    	if (!sf_getbool("inv",&inv)) inv=adj; 
	/* if y, perform inverse operation */
    	if (!sf_getbool ("verb",&verb)) verb=false;
	/* verbosity flag */

    	/* read input file parameters */
    	if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
	/* number of samples in time axis */
    	if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
	/* interval of time axis */
    	if (!sf_histfloat(in,"o1",&t0)) t0=0.;
	/* origin of time axis */
	nfft=2*kiss_fft_next_fast_size(nt);

    	if (adj) { // m(tau,p)=sum_{i=0}^{nx} d(t=tau+p*x_i,x_i)
		if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
		/* number of offset if the input in the data domain */

		/* specify slope axis */
		if (!sf_getint  ("np",&np)) sf_error("Need np=");
		/* number of p values (if adj=y) */
		if (!sf_getfloat("dp",&dp)) sf_error("Need dp=");
		/* p sampling (if adj=y) */
		if (!sf_getfloat("p0",&p0)) sf_error("Need p0=");
		/* p origin (if adj=y) */

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
	cdd=sf_complexalloc2(nfft, nx);
	cmm=sf_complexalloc2(nfft, np);
    	tmp=(fftwf_complex*)fftwf_malloc(nfft*sizeof(fftwf_complex));
    	fft1=fftwf_plan_dft_1d(nfft,tmp,tmp,FFTW_FORWARD,FFTW_MEASURE);	
   	ifft1=fftwf_plan_dft_1d(nfft,tmp,tmp,FFTW_BACKWARD,FFTW_MEASURE);

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
	for(ix=0; ix<nx; ix++)
	for(it=0; it<nfft; it++)
		cdd[ix][it]=sf_cmplx(0,0);
	for(ip=0; ip<np; ip++)
	for(it=0; it<nfft; it++)
		cmm[ip][it]=sf_cmplx(0,0);

    	if (!sf_getbool("parab",&par)) par=false;
	/* if y, parabolic Radon transform */
    	if (!sf_getfloat("x0",&x0)) x0=1.;   
	/* reference offset */

	for (ix=0; ix < nx; ix++)/* normalize offsets */
	{
		if (par) xx[ix] *= xx[ix]/(x0*x0);
		else if (x0!=1.) xx[ix] /= x0;
	}

	if(adj){// m(tau,p)=sum_{i=0}^{nx} d(t=tau+p*x_i,x_i)
		for(ix=0; ix<nx; ix++)
		for(it=0; it<nt; it++)
			cdd[ix][it]=sf_cmplx(dd[ix][it],0);

		for(ix=0; ix<nx; ix++) // loop over offsets
		{
			for(it=0; it<nfft; it++) tmp[it] = cdd[ix][it];
		 	fftwf_execute(fft1);// FFT: dd-->cdd
			for(it=0; it<nfft; it++) cdd[ix][it]=tmp[it];
		}

		for(ip=0; ip<np; ip++) // loop over slopes
		{
			for(it=0; it<=(nfft-1)/2; it++) 
			{
				tmpc=sf_cmplx(0,0);
				for(ix=0; ix<nx; ix++) 
					tmpc+=cexpf(sf_cmplx(0,2.*SF_PI*it*p[ip]*xx[ix]/nx))*cdd[ix][it];
				cmm[ip][it]=tmpc;
			}
			for(it=(nfft+1)/2; it<nfft; it++)//X(N-k)=X^*(k)
				cmm[ip][it]=conjf(cmm[ip][nfft-it]);

			for(it=0; it<nfft; it++) tmp[it] = cmm[ip][it];
		 	fftwf_execute(ifft1); // IFFT: cmm-->mm
			for(it=0; it<nt; it++) mm[ip][it]=crealf(tmp[it]);
		}

		sf_floatwrite(mm[0], nt*np, out);
	}else{// d(t,h)=sum_{i=0}^{np} m(tau=t-p_i*h,p_i)
		for(ip=0; ip<np; ip++)
		for(it=0; it<nt; it++)
			cmm[ip][it]=sf_cmplx(mm[ip][it],0);

		for(ip=0; ip<np; ip++) // loop over slopes
		{
			for(it=0; it<nfft; it++) tmp[it]=cmm[ip][it];
		 	fftwf_execute(fft1);// FFT: mm-->cmm
			for(it=0; it<nfft; it++) cmm[ip][it]=tmp[it];
		}

		for(ix=0; ix<nx; ix++) // loop over offsets
		{
			for(it=0; it<=(nfft-1)/2; it++) 
			{
				tmpc=sf_cmplx(0,0);
				for(ip=0; ip<np; ip++)
					tmpc+=cexpf(sf_cmplx(0,-2.*SF_PI*it*p[ip]*xx[ix]/nx))*cmm[ip][it];
				cdd[ix][it]=tmpc;
			}
			for(it=(nfft+1)/2; it<nfft; it++)//X(N-k)=X^*(k)
				cdd[ix][it]=conjf(cdd[ix][nfft-it]);

			for(it=0; it<nfft; it++) tmp[it]=cdd[ix][it];
		 	fftwf_execute(ifft1);// IFFT: cmm-->mm
			for(it=0; it<nt; it++) dd[ix][it]=crealf(tmp[it]);
		}

		sf_floatwrite(dd[0], nt*nx, out);
	}

	free(p);
	free(xx);
	free(*dd); free(dd);
	free(*mm); free(mm);
	free(*cdd); free(cdd);
	free(*cmm); free(cmm);
	fftwf_free(tmp);
	fftwf_destroy_plan(fft1);
    	fftwf_destroy_plan(ifft1);

    	exit(0);
}

