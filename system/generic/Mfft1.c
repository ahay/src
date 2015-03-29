/* Fast Fourier Transform along the first axis. */
/*
  Copyright (C) 2004 University of Texas at Austin
  Copyright (C) 2014 Colorado School of Mines

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
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef SF_HAS_FFTW
#include <fftw3.h>
#endif

int main (int argc, char *argv[])
{
    bool  verb,inv, sym, opt;
    int   n1, nw, nt,i1;
    float o1;
    float d1, dw;
    float shift;
    float wght;
    char *label;
    sf_file Fin=NULL, Fou=NULL;

#ifdef SF_HAS_FFTW
#ifdef SF_HAS_FFTW_OMP
    fftwf_plan    ompcfg;
#else
    fftwf_plan    *ompcfg;
#endif
#else
    kiss_fftr_cfg *ompcfg;
#endif

    long int nbuf, ibuf, left, n2buf;
    int            ompnth=1; /* number of threads */
    int            ompith=0; /* thread index */
    float        **ompP;
    kiss_fft_cpx **ompQ;
    kiss_fft_cpx  *ompE;

    float  memsize; /* in Mb */

    /*------------------------------------------------------------*/
    sf_init(argc, argv);
    ompnth=omp_init();

    Fin  = sf_input ("in");
    Fou = sf_output("out");

    if(! sf_getbool("verb",&verb)) verb=false; /* Verbosity flag */
    if (!sf_getbool("inv",&inv))    inv=false; /* y, perform inverse transform */
    if (!sf_getbool("sym",&sym))    sym=false; /* y, symmetric scaling for Hermitian FFT */
    if (!sf_getbool("opt",&opt))    opt=true;  /* y, determine optimal size for efficiency */

    if (inv) {
	if (SF_COMPLEX != sf_gettype(Fin)) sf_error("Need complex input");
	sf_settype (Fou,SF_FLOAT);
    } else {
	if (SF_FLOAT   != sf_gettype(Fin)) sf_error("Need float input");
	sf_settype (Fou,SF_COMPLEX);
    }

    if (!inv) {
	if (!sf_histint  (Fin,"n1",&n1)) n1=1;
	if (!sf_histfloat(Fin,"d1",&d1)) d1=1.;
	if (!sf_histfloat(Fin,"o1",&o1)) o1=0.;

	/* determine wavenumber sampling (for real to complex FFT) */
	nt = opt? 2*kiss_fft_next_fast_size((n1+1)/2): n1;
	if (nt%2) nt++;
	nw = nt/2+1;
	dw = 1./(nt*d1);

	sf_putint  (Fou,"n1",nw);
	sf_putfloat(Fou,"o1",0.);
	sf_putfloat(Fou,"d1",dw);
	sf_putfloat(Fou,"fft_o1",o1);
	sf_putfloat(Fou,"fft_n1",n1);

	/* fix label */
	if (NULL != (label = sf_histstring(Fin,"label1"))) {
	    sf_putstring(Fou,"fft_label1",label);
	    if (!sf_fft_label(1,label,Fou))
		sf_putstring(Fou,"label1","Wavenumber");
	}
    } else {
	if (!sf_histint  (Fin,"n1",&nw)) sf_error("No n1= in input");
	if (!sf_histfloat(Fin,"d1",&dw)) sf_error("No d1= in input");
	if (!sf_histfloat(Fin,"fft_o1",&o1)) o1=0.; 

	nt = 2*(nw-1);
	d1 = 1./(nt*dw);

	if (!opt || !sf_histint(Fin,"fft_n1",&n1)) n1 = nt;

	sf_putint  (Fou,"n1",n1);
	sf_putfloat(Fou,"d1",d1);
	sf_putfloat(Fou,"o1",o1);

	/* fix label */
	if (NULL != (label = sf_histstring(Fin,"fft_label1"))) {
	    sf_putstring(Fou,"label1",label);
	} else if (NULL != (label = sf_histstring(Fin,"label1"))) {
	    (void) sf_fft_label(1,label,Fou);
	}
    }
    sf_fft_unit(1,sf_histstring(Fin,"unit1"),Fou);

    /*------------------------------------------------------------*/
    /* Hermitian weight */
    wght = sym ? 1.0f/sqrtf((float) nt): 1.0f/nt;

    /*------------------------------------------------------------*/
    /* usable memory (Mb) */
    if(!sf_getfloat("memsize",&memsize)) memsize=1000.0;
    if(verb) sf_warning("memsize=%g",memsize);

    n2buf=1;
    while(n2buf/1024.*n1/1024.*SF_FLOAT < memsize) n2buf++;
    n2buf=SF_MIN(n2buf,sf_leftsize(Fin,1));
    /* sf_warning("n2buf=%ld",n2buf); */

    /*------------------------------------------------------------*/
    if(verb) sf_warning("allocate arrays %d %d",n2buf,ompnth);

    ompP =                  sf_floatalloc2  (nt,n2buf);
    ompQ = (kiss_fft_cpx**) sf_complexalloc2(nw,n2buf);
    ompE = (kiss_fft_cpx* ) sf_complexalloc (ompnth);

    /*------------------------------------------------------------*/
    if(verb) sf_warning("init FFT");

#ifdef SF_HAS_FFTW
#ifdef SF_HAS_FFTW_OMP
    fftwf_init_threads();
    fftwf_plan_with_nthreads(ompnth);
    if (inv) {
	ompcfg = fftwf_plan_many_dft_c2r(1, &nt, n2buf,
					 (fftwf_complex *) ompQ[0], NULL, 1, nw,
					 ompP[0], NULL, 1, nt,
					 FFTW_ESTIMATE);
    } else {
	ompcfg = fftwf_plan_many_dft_r2c(1, &nt, n2buf,
					 ompP[0], NULL, 1, nt,
					 (fftwf_complex *) ompQ[0], NULL, 1, nw,
					 FFTW_ESTIMATE);
    }
    if (NULL == ompcfg) sf_error("FFTW failure.");
#else
    ompcfg  = (fftwf_plan*) sf_alloc(n2buf,sizeof(fftwf_plan));
    for(ibuf=0; ibuf<n2buf; ibuf++)
	if (inv) ompcfg[ibuf] = fftwf_plan_dft_c2r_1d(nt,             (fftwf_complex *) ompQ[ibuf], ompP[ibuf],FFTW_ESTIMATE);
	else     ompcfg[ibuf] = fftwf_plan_dft_r2c_1d(nt, ompP[ibuf], (fftwf_complex *) ompQ[ibuf],            FFTW_ESTIMATE);    
#endif
#else
    ompcfg  = (kiss_fftr_cfg*) sf_alloc(ompnth,sizeof(kiss_fftr_cfg));
    for(ompith=0; ompith<ompnth; ompith++)
	ompcfg[ompith] = kiss_fftr_alloc(nt,inv?1:0,NULL,NULL);
#endif

    /*------------------------------------------------------------*/
    ompith=0;
    nbuf = n2buf;
    for (left=sf_leftsize(Fin,1); left>=0; left -= nbuf) {
	if(verb) sf_warning("%ld %ld;",left,nbuf);

	/* buffer size */
	nbuf=SF_MIN(left,nbuf);

	if (!inv) { /* FORWARD TRANSFORM */
	    for(ibuf=0; ibuf<nbuf; ibuf++)
		sf_floatread (ompP[ibuf],n1,Fin);

#ifdef SF_HAS_FFTW_OMP
	    for(ibuf=0; ibuf<nbuf; ibuf++) {
		if (sym) for(i1=0;  i1<n1; i1++) ompP[ibuf][i1] *= wght;
		;        for(i1=n1; i1<nt; i1++) ompP[ibuf][i1]  = 0.0;
	    }
	    fftwf_execute(ompcfg);
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(static)			\
    private(ibuf,i1,ompith,shift)				\
    shared( nbuf,n1,nt,nw,dw,o1,wght,ompP,ompQ,ompE,ompcfg)
#endif
	    for(ibuf=0; ibuf<nbuf; ibuf++) {
#ifdef _OPENMP
		ompith = omp_get_thread_num();
#endif

#ifndef SF_HAS_FFTW_OMP
		if (sym) for(i1=0;  i1<n1; i1++) ompP[ibuf][i1] *= wght;
		;        for(i1=n1; i1<nt; i1++) ompP[ibuf][i1]  = 0.0;
#ifdef SF_HAS_FFTW
		fftwf_execute(ompcfg[ibuf]);
#else
		kiss_fftr(ompcfg[ompith],ompP[ibuf],ompQ[ibuf]);
#endif
#endif

		if (0. != o1) { shift = -2.0*SF_PI*dw*o1;
		    for (i1=0; i1<nw; i1++) {
			ompE[ompith].r = cosf(shift*i1);
			ompE[ompith].i = sinf(shift*i1);
			ompQ[ibuf][i1]=sf_cmul(ompQ[ibuf][i1],ompE[ompith]);
		    }
		}

	    }

	    for(ibuf=0; ibuf<nbuf; ibuf++)
		sf_floatwrite((float*) ompQ[ibuf],2*nw,Fou);
	} else {
	    for(ibuf=0; ibuf<nbuf; ibuf++)
		sf_floatread( (float*) ompQ[ibuf],2*nw,Fin);

#ifdef SF_HAS_FFTW_OMP
	    for(ibuf=0; ibuf<nbuf; ibuf++) {
		if (0. != o1) { shift = +2.0*SF_PI*dw*o1;
		    for(i1=0; i1<nw; i1++) {
			ompE[ompith].r = cosf(shift*i1);
			ompE[ompith].i = sinf(shift*i1);
			ompQ[ibuf][i1]=sf_cmul(ompQ[ibuf][i1],ompE[ompith]);
		    }
		}
	    }
	    fftwf_execute(ompcfg);
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(static)			\
    private(ibuf,i1,ompith,shift)				\
    shared( nbuf,n1,nw,dw,o1,wght,ompP,ompQ,ompE,ompcfg)
#endif
	    for(ibuf=0; ibuf<nbuf; ibuf++) {
#ifdef _OPENMP
		ompith = omp_get_thread_num();
#endif

#ifndef SF_HAS_FFTW_OMP
		if (0. != o1) { shift = +2.0*SF_PI*dw*o1;
		    for(i1=0; i1<nw; i1++) {
			ompE[ompith].r = cosf(shift*i1);
			ompE[ompith].i = sinf(shift*i1);
			ompQ[ibuf][i1]=sf_cmul(ompQ[ibuf][i1],ompE[ompith]);
		    }
		}
#ifdef SF_HAS_FFTW
		fftwf_execute(ompcfg[ibuf]);
#else
		kiss_fftri(ompcfg[ompith],ompQ[ibuf],ompP[ibuf]);
#endif
#endif

		for(i1=0; i1<n1; i1++) ompP[ibuf][i1] *= wght;
	    }

	    for(ibuf=0; ibuf<nbuf; ibuf++)
		sf_floatwrite (ompP[ibuf],n1,Fou);
	}
    }
    if(verb) sf_warning(".");

    /*------------------------------------------------------------*/

    if(verb) sf_warning("deallocate arrays %d %d",n2buf,ompnth);

    free(*ompP); free(ompP);
    free(*ompQ); free(ompQ);
    ;            free(ompE);

    exit (0);
}
