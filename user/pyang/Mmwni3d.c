/* bandlimited minimum weighted-norm interpolation (MWNI) 
 implemented with conjugate gradient least squares (CGLS) method
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

void lowpass_filter(fftwf_complex *tmp, int n1, int n2, int n3)
/*< only low-wavenumber preserved >*/
{
	for(int i3=0;i3<n3;i3++){
		for(int i2=0;i2<n2;i2++){			
	    	    if( (i3>SF_NINT(n3*0.15) && i3<SF_NINT(n3*0.85)) ||
			(i2>SF_NINT(n2*0.15) && i2<SF_NINT(n2*0.85)) ) {	
			for(int i1=0; i1<n1; i1++) 
				tmp[i1+n1*(i2+n2*i3)]=0.0;
		    }	
		}
	}
}


int main(int argc, char* argv[])
{
    /* input and output variables */
    bool verb;
    int niter; 
    float tol;
    sf_file Fin, Fout, Fmask;/* mask and I/O files*/ 

    /* define temporary variables */
    int n1,n2,n3;
    float *din, *mask;

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
    float g0, gnp, gn, beta, alpha;
    fftwf_complex *gm;
    sf_complex *m, *r, *gr, *s;
    fftwf_plan fft3, ifft3;/* execute plan for FFT and IFFT */

    m=sf_complexalloc(n1*n2*n3);
    r=sf_complexalloc(n1*n2*n3);
    gr=sf_complexalloc(n1*n2*n3);
    s=sf_complexalloc(n1*n2*n3);
    gm=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*n1*n2*n3);
    fft3=fftwf_plan_dft_3d(n1,n2,n3,gm,gm,FFTW_FORWARD,FFTW_MEASURE);	
    ifft3=fftwf_plan_dft_3d(n1,n2,n3,gm,gm,FFTW_BACKWARD,FFTW_MEASURE);

    // initialization
    for(int i=0; i<n1*n2*n3; i++) 
    {
	r[i]=-din[i];
	m[i]=0.0;
	gm[i]=0.0;
	gr[i]=0.0;
	s[i]=0.0;
    }
    g0=gnp=gn=beta=alpha=0;

    // conjugate gradient loop
    for(int iter=0;iter<niter;iter++)
    {
	apply_mask(mask, r, gm, n1, n2, n3); 		// gm=M.*r
	fftwf_execute(ifft3);				// gm=ifft3(M.*r);
	for(int i=0;i<n1*n2*n3;i++) gm[i]/=sqrtf(n1*n2*n3);
	lowpass_filter(gm, n1, n2, n3);		// gm=p.*fft3(M.*r);
	gn=cblas_scnrm2(n1*n2*n3, gm, 1);		// gn=sum(abs(gm(:)).^2);

	if (iter == 0) {
	    beta=0.0;
	    g0=gn;
	} else {
	    beta=gn/gnp;
	    if(beta<tol ||gn/g0<tol) break;
	}
	gnp=gn;

	for(int i=0;i<n1*n2*n3;i++){
	    s[i]=gm[i]+beta*s[i]; 			// s=gm+beta*s
	    gm[i]=s[i];					// copy s(:) to gm(:)
	}

	fftwf_execute(fft3);				// gm=fft3(s);
	for(int i=0;i<n1*n2*n3;i++) gm[i]/=sqrtf(n1*n2*n3);
	apply_mask(mask, gm, gr, n1, n2, n3); 		// gr=M.*fft3(gm);

    	alpha=-gn/cblas_scnrm2 (n1*n2*n3, gr, 1); 	// alpha=-gn/sum(abs(gr(:)).^2);
	for(int i=0;i<n1*n2*n3;i++){ 
	    m[i]+=alpha*s[i];				// m=m+alpha*s;
	    r[i]+=alpha*gr[i];				// r=r+alpha*gr;
	}	
	if (verb) sf_warning("iteration %d of %d res %g", iter+1, niter, cblas_scnrm2 (n1*n2*n3, r, 1));
    }
    for(int i=0; i<n1*n2*n3; i++) gm[i]=m[i];    	// copy m(:) to g(:)
    fftwf_execute(fft3);			 	// gm=fft3(m);
    for(int i=0; i<n1*n2*n3; i++) din[i]=crealf(gm[i])/sqrtf(n1*n2*n3);	// drec=real(ifft2(m)*n1*n2);

    free(m);
    free(r);
    free(gr);
    free(s);
    fftwf_free(gm);
    fftwf_destroy_plan(fft3);
    fftwf_destroy_plan(ifft3);
    /* ============================================================================= */

    sf_floatwrite(din,n1*n2*n3,Fout); /* output reconstructed seismograms */

    sf_close();
    exit(0);
}
