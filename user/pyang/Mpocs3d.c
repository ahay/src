/* POCS for 3D missing data interpolation
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
#include <complex.h>
#include <fftw3.h>

#ifdef _OPENMP
#include <omp.h>
#endif


int main(int argc, char* argv[])
{
    /* input and output variables */
    bool verb;
    int niter; 
    sf_file Fin=NULL,Fout=NULL, Fmask=NULL;/* mask and I/O files*/ 

    /* define temporary variables */
    int n1,n2,n3;
    float *din=NULL, *mask=NULL, *dout=NULL;

    sf_init(argc,argv);	/* Madagascar initialization */
    omp_init(); 	/* initialize OpenMP support */

    /* setup I/O files */
    Fin=sf_input("in");	/* read the data to be interpolated */
    Fmask=sf_input("mask");  	/* read the 2-D mask for 3-D data */
    Fout=sf_output("out"); 	/* output the reconstructed data */
 
    if(!sf_getbool("verb",&verb)) 	verb=false;
    /* verbosity */
    if (!sf_getint("niter",&niter)) 	niter=100;
    /* total number of POCS iterations */

    /* Read the data size */
    if (!sf_histint(Fin,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(Fin,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(Fin,"n3",&n3)) sf_error("No n3= in input");

    /* allocate data and mask arrays */
    din=sf_floatalloc(n1*n2*n3); sf_floatread(din,n1*n2*n3,Fin);
    dout=sf_floatalloc(n1*n2*n3);
    if (NULL != sf_getstring("mask")){
	mask=sf_floatalloc(n2*n3);
	sf_floatread(mask,n2*n3,Fmask);
    }

    /********************* 3-D POCS interpolation *********************/
    fftwf_plan p1,p2;/* execute plan for FFT and IFFT */
    int i1,i2,i3, iter;
    float thr;		
    fftwf_complex *dprev,*dcurr;

    dprev=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*n1*n2*n3);
    dcurr=(fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex)*n1*n2*n3);
    p1=fftwf_plan_dft_3d(n1,n2,n3,dprev,dcurr,FFTW_FORWARD,FFTW_MEASURE);	
    p2=fftwf_plan_dft_3d(n1,n2,n3,dcurr,dprev,FFTW_BACKWARD,FFTW_MEASURE);

    for(i1=0; i1<n1*n2*n3; i1++) dprev[i1]=din[i1];

    /* POCS iterations */
    for(iter=1; iter<=niter; iter++)
    {
	fftwf_execute(p1);/* FFT */

	/* find the thresholds sequence */
	float mmax=cabs(dcurr[0]);
	for(i1=1; i1<n1*n2*n3; i1++){
	    mmax=(cabs(dcurr[i1])>mmax)?cabs(dcurr[i1]):mmax;
	}
	thr=powf(0.005,(iter-1.0)/(niter-1.0))*mmax;

	/* perform hard thresholding */
#ifdef _OPENMP
#pragma omp parallel for	\
	private(i1)		\
	shared(dcurr,n1,n2,n3,thr)
#endif
	for(i1=0;i1<n1*n2*n3;i1++) {
	    dcurr[i1]=dcurr[i1]*(cabs(dcurr[i1])>thr?1.:0.);
	}

	fftwf_execute(p2);/* unnormalized IFFT */

#ifdef _OPENMP
#pragma omp parallel for	\
	private(i1)		\
	shared(dprev,n1,n2,n3)
#endif
	/* adjointness needs scaling with factor 1.0/(n1*n2*n3) */	
	for(i1=0; i1<n1*n2*n3;i1++) dprev[i1]=dprev[i1]/(n1*n2*n3);
	
	/* update d_rec: d_rec = d_obs+(1-M)*A T{ At(d_rec) } */
	for(i3=0;i3<n3;i3++)	
	    for(i2=0;i2<n2;i2++){
		if (mask[i2+i3*n2]!=0.){			
		    for(i1=0; i1<n1; i1++)
			dprev[i1+n1*(i2+n2*i3)]=din[i1+n1*(i2+n2*i3)];
		    }
		}
    }

    /* take the real part */
    for(i1=0;i1<n1*n2*n3; i1++) dout[i1]=creal(dprev[i1]);
	
    fftwf_destroy_plan(p1);
    fftwf_destroy_plan(p2);
    fftwf_free(dprev);
    fftwf_free(dcurr);

    /* write reconstructed seismic data to output */
    sf_floatwrite(dout,n1*n2*n3,Fout);

    sf_close();
    exit(0);
}
