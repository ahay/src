/* Two-step POCS interpolation using a general Lp-norm optimization
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
*/
#include <rsf.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif


void apply_mask(float *mask, sf_complex *before, sf_complex *after, int n1, int n2, int n3)
/*< apply the mask >*/
{
	for(int i3=0;i3<n3;i3++)	
	for(int i2=0;i2<n2;i2++)
	{
		if (mask[i2+i3*n2]){			
		    for(int i1=0; i1<n1; i1++) after[i1+n1*(i2+n2*i3)]=before[i1+n1*(i2+n2*i3)];
		}else{		
		    for(int i1=0; i1<n1; i1++) after[i1+n1*(i2+n2*i3)]=0.0;			
		}
	}
}

void lowpass_filter(sf_complex *tmp, int n1, int n2, int n3)
/*< only low-wavenumber preserved >*/
{
	for(int i3=0;i3<n3;i3++){
	    if (i3>SF_NINT(n3*0.2) && i3<SF_NINT(n3*0.8)){	
		for(int i2=0;i2<n2;i2++){			
	    	    if (i2>SF_NINT(n2*0.2) && i2<SF_NINT(n2*0.8)){	
			for(int i1=0; i1<n1; i1++) 
				tmp[i1+n1*(i2+n2*i3)]=0.0;
		    }	
		}
     	    }
	}
}

float compute_squaredsum(sf_complex *x, int n)
/*< sum |x_i|^2>*/
{
	float sum=0;
	for(int i=0;i<n;i++) sum+=cabsf(x[i])*cabsf(x[i]);
	return sum;
}

void compute_xpay(sf_complex *y, sf_complex *x, float a, int n)
/*< x plus a y: y(:)=x(:)+a*y(:) >*/
{
	int i;
#ifdef _OPENMP
#pragma omp parallel for	\
	schedule(dynamic) 	\
	private(i)		\
	shared(a)
#endif
	for(i=0;i<n;i++) y[i]=x[i]+a*y[i];
}


void compute_axpy(sf_complex *y, sf_complex *x, float a, int n)
/*< a x plus y: y(:)=y(:)+a*x(:) >*/
{
	int i;
#ifdef _OPENMP
#pragma omp parallel for	\
	schedule(dynamic) 	\
	private(i)		\
	shared(a)
#endif
	for(i=0;i<n;i++) y[i]+=a*x[i];
}

int main(int argc, char* argv[])
{
    /* input and output variables */
    bool verb;
    int niter; 
    float tol;
    sf_file Fin=NULL,Fout=NULL, Fmask=NULL;/* mask and I/O files*/ 

    /* define temporary variables */
    int n1,n2,n3;
    float *din=NULL, *mask=NULL;

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

    /* Read the data size */
    if (!sf_histint(Fin,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(Fin,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(Fin,"n3",&n3)) sf_error("No n3= in input");

    /* allocate data and mask arrays */
    din=sf_floatalloc(n1*n2*n3); sf_floatread(din,n1*n2*n3,Fin);
    if (NULL != sf_getstring("mask")){
	mask=sf_floatalloc(n2*n3);
	sf_floatread(mask,n2*n3,Fmask);
    }

    /*=========================== conjugate gradient iterations ================= */
    float rho0, rho_old, rho_new, beta, alpha;
    fftwf_complex *m, *r, *sm, *sr, *gm, *gr;
    fftwf_plan fft3,ifft3;/* execute plan for FFT and IFFT */

    m=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*n1*n2*n3);
    r=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*n1*n2*n3);
    gm=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*n1*n2*n3);
    gr=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*n1*n2*n3);
    sm=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*n1*n2*n3);
    sr=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*n1*n2*n3);
    fft3=fftwf_plan_dft_3d(n1,n2,n3,gm,gm,FFTW_FORWARD,FFTW_MEASURE);	
    ifft3=fftwf_plan_dft_3d(n1,n2,n3,gm,gm,FFTW_BACKWARD,FFTW_MEASURE);

    // initialization
    for(int i=0; i<n1*n2*n3; i++) 
    {
	r[i]=-din[i];
	m[i]=0.0;
	sm[i]=0.0;
	sr[i]=0.0;
	gm[i]=0.0;
	gr[i]=0.0;
    }
    // conjugate gradient loop
    for(int iter=0;iter<niter;iter++)
    {
	apply_mask(mask, r, gm, n1, n2, n3); 		// gm=M.*r
	fftwf_execute(fft3);				// gm=fft3(M.*r);
	lowpass_filter(gm, n1, n2, n3);			// gm=p.*fft3(M.*r);
	fftwf_execute(ifft3);				// gm=ifft3(gm)*n1*n2;
	apply_mask(mask, gm, gr, n1, n2, n3); 		// gr=M.*ifft3(gm)*n1*n2;
	rho_new=compute_squaredsum(gm,n1*n2*n3);	// rho1=sum(abs(gm(:)).^2);

	if (iter == 0) {
	    beta=0.0;
	    rho0=rho_new;
	} else {
	    beta=rho_new/rho_old;
	    if(beta<tol || rho_new/rho0<tol) break;

   	    compute_xpay(sm, gm, beta, n1*n2*n3);	// sm=gm+beta*sm;
   	    compute_xpay(sr, gr, beta, n1*n2*n3);	// sr=gr+beta*sr;
	    float sr2=compute_squaredsum(sr,n1*n2*n3);
    	    alpha=-rho_new/sr2; 		 	//alpha=-rho1/sum(abs(sr(:)).^2);
   	    compute_axpy(m, sm, alpha, n1*n2*n3);	// m=m+alpha*sm;
   	    compute_axpy(r, sr, alpha, n1*n2*n3);	// r=r+alpha*sr;
	}
	rho_old=rho_new;
    }
    for(int i=0; i<n1*n2*n3; i++) gm[i]=m[i];    	// copy m(:) to g(:)
    fftwf_execute(ifft3);			 	// gm=ifft3(m)*n1*n2;
    for(int i=0; i<n1*n2*n3; i++) din[i]=crealf(gm[i]);	// drec=real(ifft2(m)*n1*n2);

    fftwf_free(m);
    fftwf_free(r);
    fftwf_free(gm);
    fftwf_free(gr);
    fftwf_free(sm);
    fftwf_free(sr);
    fftwf_destroy_plan(fft3);
    fftwf_destroy_plan(ifft3);

    /* write reconstructed seismic data to output */
    sf_floatwrite(din,n1*n2*n3,Fout);

    sf_close();
    exit(0);
}
