/* 1D Fourier transform */
/*
  Copyright (C) 2016 Yangkang Chen
*/

#include<rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef SF_HAS_FFTW
#include <fftw3.h>
#endif

#include"fft1.h"

	 	    
void fft1(float **ompP /*input data*/, 
		kiss_fft_cpx **ompQ,
		int n2,
		int n1,
		float d1,
		float o1, 		
		int nt,
		int nw,
		float dw,
		bool sym, 
		bool opt,
		bool verb, 		
		bool inv)	
/*<1D forward and inverse Fourier transform>*/    
{

    float shift;
    float wght;
    int   ompnth=1; /* number of threads */
    int   ompith=0; /* thread index */
    int i1,i2;
    kiss_fft_cpx  *ompE;
            
    ompnth=omp_init();
    
    
#ifdef SF_HAS_FFTW
#ifdef SF_HAS_FFTW_OMP
    fftwf_plan    ompcfg;
#else
    fftwf_plan    *ompcfg;
#endif
#else
    kiss_fftr_cfg *ompcfg;
#endif

    /*------------------------------------------------------------*/
    /* Hermitian weight */
    wght = sym ? 1.0f/sqrtf((float) nt): 1.0f/nt;


    /*------------------------------------------------------------*/
    if(verb) sf_warning("allocate arrays %d %d",n2,ompnth);


    ompE = (kiss_fft_cpx* ) sf_complexalloc (ompnth);

    /*------------------------------------------------------------*/
    if(verb) sf_warning("init FFT");

#ifdef SF_HAS_FFTW
#ifdef SF_HAS_FFTW_OMP
    fftwf_init_threads();
    fftwf_plan_with_nthreads(ompnth);
    if (inv) {
	ompcfg = fftwf_plan_many_dft_c2r(1, &nt, n2,
					 (fftwf_complex *) ompQ[0], NULL, 1, nw,
					 ompP[0], NULL, 1, nt,
					 FFTW_ESTIMATE);
    } else {
	ompcfg = fftwf_plan_many_dft_r2c(1, &nt, n2,
					 ompP[0], NULL, 1, nt,
					 (fftwf_complex *) ompQ[0], NULL, 1, nw,
					 FFTW_ESTIMATE);
    }
    if (NULL == ompcfg) sf_error("FFTW failure.");
#else
    ompcfg  = (fftwf_plan*) sf_alloc(n2,sizeof(fftwf_plan));
    for(i2=0; i2<n2; i2++)
	if (inv) ompcfg[i2] = fftwf_plan_dft_c2r_1d(nt,             (fftwf_complex *) ompQ[i2], ompP[i2],FFTW_ESTIMATE);
	else     ompcfg[i2] = fftwf_plan_dft_r2c_1d(nt, ompP[i2], (fftwf_complex *) ompQ[i2],            FFTW_ESTIMATE);    
#endif
#else
    ompcfg  = (kiss_fftr_cfg*) sf_alloc(ompnth,sizeof(kiss_fftr_cfg));
    for(ompith=0; ompith<ompnth; ompith++)
	ompcfg[ompith] = kiss_fftr_alloc(nt,inv?1:0,NULL,NULL);
    ompith=0;
#endif

	if (!inv) { /* FORWARD TRANSFORM */


#ifdef SF_HAS_FFTW_OMP
	    for(i2=0; i2<n2; i2++) {
		if (sym) for(i1=0;  i1<n1; i1++) ompP[i2][i1] *= wght;
		;        for(i1=n1; i1<nt; i1++) ompP[i2][i1]  = 0.0;
	    }
	    fftwf_execute(ompcfg);
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(static)			\
    private(i2,i1,ompith,shift)				\
    shared( n2,n1,nt,nw,dw,o1,wght,ompP,ompQ,ompE,ompcfg)
#endif
	    for(i2=0; i2<n2; i2++) {
#ifdef _OPENMP
		ompith = omp_get_thread_num();
#endif

#ifndef SF_HAS_FFTW_OMP
		if (sym) for(i1=0;  i1<n1; i1++) ompP[i2][i1] *= wght;
		;        for(i1=n1; i1<nt; i1++) ompP[i2][i1]  = 0.0;
#ifdef SF_HAS_FFTW
		fftwf_execute(ompcfg[i2]);
#else
		kiss_fftr(ompcfg[ompith],ompP[i2],ompQ[i2]);
#endif
#endif

		if (0. != o1) { shift = -2.0*SF_PI*dw*o1;
		    for (i1=0; i1<nw; i1++) {
			ompE[ompith].r = cosf(shift*i1);
			ompE[ompith].i = sinf(shift*i1);
			ompQ[i2][i1]=sf_cmul(ompQ[i2][i1],ompE[ompith]);
		    }
		}

	    }

	} else { /* INVERSE TRANSFORM */

#ifdef SF_HAS_FFTW_OMP
	    ompith=0;

	    for(i2=0; i2<n2; i2++) {
		if (0. != o1) { shift = +2.0*SF_PI*dw*o1;
		    for(i1=0; i1<nw; i1++) {
			ompE[ompith].r = cosf(shift*i1);
			ompE[ompith].i = sinf(shift*i1);
			ompQ[i2][i1]=sf_cmul(ompQ[i2][i1],ompE[ompith]);
		    }
		}
	    }
	    fftwf_execute(ompcfg);
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(static)			\
    private(i2,i1,ompith,shift)				\
    shared( n2,n1,nw,dw,o1,wght,ompP,ompQ,ompE,ompcfg)
#endif
	    for(i2=0; i2<n2; i2++) {
#ifdef _OPENMP
		ompith = omp_get_thread_num();
#endif

#ifndef SF_HAS_FFTW_OMP
		if (0. != o1) { shift = +2.0*SF_PI*dw*o1;
		    for(i1=0; i1<nw; i1++) {
			ompE[ompith].r = cosf(shift*i1);
			ompE[ompith].i = sinf(shift*i1);
			ompQ[i2][i1]=sf_cmul(ompQ[i2][i1],ompE[ompith]);
		    }
		}
#ifdef SF_HAS_FFTW
		fftwf_execute(ompcfg[i2]);
#else
		kiss_fftri(ompcfg[ompith],ompQ[i2],ompP[i2]);
#endif
#endif

		for(i1=0; i1<n1; i1++) ompP[i2][i1] *= wght;
	    }

	}
}

