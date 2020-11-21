/* Marchenko-Wapenaar-Broggini iterative scheme

sfmarchenko < downgoing.rsf refl=REFL_000.rsf conj=y verb=n Gtot=y niter=21 nshots=401 scale=1 eps=1e-4 shift=5 Gm=Gm.rsf G=G.rsf> Gp.rsf

======= INPUTS ============

p0plus.rsf: initial downgoing wavefield

REFL_000.rsf: Fourier transform of the reflection response

======= PARAMETERS ========

conj  = [y]/n	- complex-conjugation of the first input (corresponds to time-reversal in time)
verb = y/[n]	- verbosity flag
twin  = y/[n]	- returns the timewindow as one of the outputs (window=window.rsf)
pandq  = y/[n]	- pandq=true: returns p and q, pandq=false returns Gp and Gm
Gtot  = y/[n]	- Gtot=true returns G=Gp+Gm
Htot  = y/[n]	- Htot=true returns H=Gp-Gm
niter  = 1		- number of iterations
nshots  = 1		- number of shots in the reflection response
scale  = 1.0	- scale factor (often due to resampling)
eps  = 1e-4		- threshold for the timewindow
shift  = 5		- shift in samples for the timewindow
*/

/*
  Copyright (C) 2012 Colorado School of Mines
  
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
#include "fft1.h"

void fft1 (float *, float *, sf_file, bool, bool, bool);

int main(int argc, char* argv[])
{

	bool verb,conj,twin,pandq,Gtot,Htot;
	
	float	*pplus0, *pplus, *pplusinv, *pplustemp, *Pplus, *Pplus_trace, *pminus, *Pminus, *Refl; 
	float   *Gp, *Gm, *G = NULL, *H = NULL;
	float	*qplus, *qplustemp, *Qplus, *Qplus_trace, *qminus, *Qminus;
	float	*window, *taper, pi;
	int		*tw;

    /* I/O files */
    sf_file Fplus;
    sf_file FRefl;
    sf_file FGp;
    sf_file FGm;
    sf_file FG = NULL;
    sf_file FH = NULL;
    sf_file Ftwin = NULL;
    sf_file Fp = NULL;
    sf_file Fq = NULL;

	char *filename1, filename2[256], filename3[256];
	
	/* Cube axes */
    sf_axis at,af,ax;

    int     nt,nf,ntr,mode,nshots,niter,len;
    int     i,it,ix,ishot,iter,i0;
	int		twc, twa, shift, n[2], rect[2], s[2];
    float   scale,eps,a,b,c,d,e,f;

	sf_triangle tr;

    /*------------------------------------------------------------*/
    /* Initialize RSF parameters 								  */
    /*------------------------------------------------------------*/
    sf_init(argc,argv);	
	
    /*------------------------------------------------------------*/
    /* Initialize OMP parameters */
    /*------------------------------------------------------------*/
	#ifdef _OPENMP
    omp_init();
	#endif	

	/*------------------------------------------------------------*/
	/* Flags 													  */
	/*------------------------------------------------------------*/
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("conj",&conj)) conj=false; /* complex conjugation (time-reversal) flag */
    if(! sf_getbool("twin",&twin)) twin=false; /* returns the timewindow as one of the outputs */
    if(! sf_getbool("pandq",&pandq)) pandq=false; /* pandq=true: returns p and q */
    if(! sf_getbool("Gtot",&Gtot)) Gtot=false; /* Gtot=true: returns G=Gp+Gm */
    if(! sf_getbool("Htot",&Htot)) Htot=false; /* Htot=true: returns H=Gp-Gm */
    if(! sf_getint("niter",&niter)) niter=1; /* number of iterations */
    if(! sf_getint("nshots",&nshots)) nshots=1; /* number of shots */
    if(! sf_getfloat("scale",&scale)) scale=1.0; /* scale factor */
	if(! sf_getfloat("eps",&eps)) eps=1e-4; /* threshold for the timewindow */
	if(! sf_getint("shift",&shift)) shift=5; /* shift in samples for the timewindow */
    
	if (verb) {
		fprintf(stderr,"This program was called with \"%s\".\n",argv[0]);
		/*fprintf(stderr,"Nr: %d Nx: %d Nt:%d\n",nr,nx,nt);*/
		
		if (argc > 1) {
			for (i = 1; i<argc; i++) {
				fprintf(stderr,"argv[%d] = %s\n", i, argv[i]);
			}
		}
		else {
			fprintf(stderr,"The command had no other arguments.\n");
    	}	
	}

    /*------------------------------------------------------------*/
    /* I/O files 												  */
    /*------------------------------------------------------------*/	
	/* "in" is the transposed version of p00plus_xxxx_xxxx.rsf
	   Dimensions of p00plus_xxxx_xxxx.rsf BEFORE sftransp: n1=ntr,n2=nt
	   Dimensions of p00plus_xxxx_xxxx.rsf BEFORE sftransp: n1=nt,n2=ntr */
	Fplus = sf_input("in");
	
	/* refl is REFL_000.rsf
	   It is used to read nf, df, of
	  Dimensions are: n1=nf,n2=ntr */
	/*FRefl = (sf_file)sf_alloc(1,sizeof(sf_file));*/	
	FRefl = sf_input("refl");

	FGp = sf_output("out");
	FGm = sf_output("Gm");
	
	if (Gtot) {
		FG  = sf_output("G");
	}
	
	if (Htot) {
		FH  = sf_output("H");
	}
	
	if (pandq) {
		Fp = sf_output("p");
		Fq = sf_output("q");
	}
	
	if (twin) {
		Ftwin = sf_output("window"); /* time window */
	}
	
	/*------------------------------------------------------------*/
	/* Axes */
	/*------------------------------------------------------------*/    
	at = sf_iaxa(Fplus,1); sf_setlabel(at,"Time"); if(verb) sf_raxa(at); /* time */
	af = sf_iaxa(FRefl,1); sf_setlabel(af,"Frequency"); if(verb) sf_raxa(af); /* frequency */
	ax = sf_iaxa(Fplus,2); sf_setlabel(ax,"r"); if(verb) sf_raxa(ax); /* space */
    
	nt = sf_n(at); 
	nf = sf_n(af); 
	ntr = sf_n(ax); 
	
	if (verb) fprintf(stderr,"nt: %d nf: %d ntr:%d\n",nt,nf,ntr);

	sf_fileclose(FRefl);

    /*------------------------------------------------------------*/
    /* Setup output data and wavefield header					  */
    /*------------------------------------------------------------*/
	sf_oaxa(FGp,at,1);
	sf_oaxa(FGp,ax,2);
	sf_oaxa(FGm,at,1);
	sf_oaxa(FGm,ax,2);
	if (Gtot) {
		sf_oaxa(FG,at,1);
		sf_oaxa(FG,ax,2);
	}
	
	if (Htot) {
		sf_oaxa(FH,at,1);
		sf_oaxa(FH,ax,2);
	}
	
	if (pandq) {
		sf_oaxa(Fp,at,1);
		sf_oaxa(Fp,ax,2);
		sf_oaxa(Fq,at,1);
		sf_oaxa(Fq,ax,2);
	}

	if (twin) {
		sf_oaxa(Ftwin,at,1);
		sf_oaxa(Ftwin,ax,2);
	}

    /*------------------------------------------------------------*/
    /* Allocate arrays											  */
    /*------------------------------------------------------------*/
	
	/* Downgoing wavefields - Time */
	pplus0 = (float *)calloc(nt*ntr,sizeof(float));
	sf_floatread(pplus0,nt*ntr,Fplus);
	pplus = (float *)calloc(nt*ntr,sizeof(float));
	memcpy(pplus,pplus0,nt*ntr*sizeof(float));
	pplustemp = (float *)calloc(nt*ntr,sizeof(float));
	pplusinv = (float *)calloc(nt*ntr,sizeof(float));
	qplus = (float *)calloc(nt*ntr,sizeof(float));
	qplustemp = (float *)calloc(nt*ntr,sizeof(float));

	/* Downgoing wavefields - Frequency */
	Pplus = (float *)calloc(2*nf*ntr,sizeof(float));
	Qplus = (float *)calloc(2*nf*ntr,sizeof(float));

	/* The three flags of fft1 are: inv, sym, and opt */
	fft1(pplus0,Pplus,Fplus,0,0,1);
	memcpy(Qplus,Pplus,2*nf*ntr*sizeof(float));
	
	/* Upgoing wavefields - Time */
	pminus = (float *)calloc(nt*ntr,sizeof(float));
	qminus = (float *)calloc(nt*ntr,sizeof(float));

	/* Downgoing wavefields - Frequency */
	Pminus = (float *)calloc(2*nf*ntr,sizeof(float));
	Qminus = (float *)calloc(2*nf*ntr,sizeof(float));

	/* This is currently NOT used */
	/* Transpose of p00plus_xxxx_xxxx */
	for (ix=0; ix<ntr; ix++) {
		for (it=0; it<nt; it++) {
			pplusinv[ix*ntr+it]=pplus0[it*ntr+ix];
		}
	}

	/* Single trace (in frequency) of the downgoing wavefield */
	Pplus_trace = (float *)calloc(2*nf,sizeof(float));
	Qplus_trace = (float *)calloc(2*nf,sizeof(float));
	
	/* Output wavefields */
	Gp = (float *)calloc(nt*ntr,sizeof(float));
	Gm = (float *)calloc(nt*ntr,sizeof(float));
	if (Gtot) {
		G = (float *)calloc(nt*ntr,sizeof(float));
	}
	
	if (Htot) {
		H = (float *)calloc(nt*ntr,sizeof(float));
	}
	
	/* Time-reversal flag */
	if (conj) {
		mode = -1;
	}
	else {
		mode = +1;
	}
    
	/* Load the reflection response into the memory */
	if (verb) fprintf(stderr,"Before loading R %d\n",2*nf*ntr);
	Refl = (float *)calloc(2*nf*ntr*nshots,sizeof(float));
	
	/* Read REFL_000.rsf */
	filename1 = sf_getstring("refl");
	/* 000.rsf are 7 characters */
	len = strlen(filename1)-7;
	/* copy the filename without 000.rsf */
	strncpy(filename2,filename1,len);
	filename2[len] = '\0';
	if (verb) fprintf(stderr,"filename2 is: %s and len is: %d\n",filename2,len);
	/*if (NULL == filename1) {
		fprintf(stderr,"Cannot read header file %s",filename1);
	}*/
  
	for (ishot=0; ishot<nshots; ishot++) {
		
	  	/* write xxx.rsf in the string filename3 */
		sprintf(filename3,"%03d.rsf",ishot);
		for (i=0; i<7; i++)
			filename2[len+i] = filename3[i];
		filename2[len+7] = '\0';
		if (verb) fprintf(stderr,"Loading %s in memory\n",filename2);
	  	FRefl = sf_input(filename2);
    	sf_floatread(&Refl[ishot*2*nf*ntr],2*nf*ntr,FRefl);
		sf_fileclose(FRefl);
		/*if (verb) fprintf(stderr,"Iteration %d\n",ishot);*/
	}

	/* Build time-window */
	tw = (int *)calloc(ntr,sizeof(int));
	window = (float *)calloc(nt*ntr,sizeof(float));
	/*memset(window,0,nt*ntr*sizeof(float));*/
    /* I am not sure why I set it to this value */
	/*for (ix=0; ix<ntr; ix++) {
		tw[ix] = nt*dt+ot+0.15; 
	}*/
	
	if (verb) fprintf(stderr,"Build the time-window\n");
	for (ix=0; ix<ntr; ix++) {
		for (it=0; it<nt; it++) {
			if (pplus0[it+ix*nt]>eps) {
				/*tw[ix] = it*dt+ot;*/
				tw[ix] = it;
				/*if (verb) fprintf(stderr,"%d %d\n",ix,it);*/
				break;
			}
		}
	}
	for (ix=0; ix<ntr; ix++) {
		twc = (int)(tw[ix]-shift);
		twa = (int)(-twc+shift+nt);
		/*if (verb) fprintf(stderr,"%d %d\n",twc,twa);*/
		for (it=0; it<nt; it++) {
			if ((it<twc) || (it>twa)) {
				window[it+ix*nt] = 1.0;
			}
		}
	}

	/* Smoothing of the window */
	/* Should I implement flags for rect and iter? */
	/* Look at Msmooth.c to understand below */
	n[0] = nt;
	n[1] = ntr;
	s[0] = 1;
	s[1] = nt;
	rect[0] = 5;
	rect[1] = 5;

	for (ix=0; ix <= 1; ix++) {
		if (rect[ix] <= 1) continue;
		tr = sf_triangle_init (rect[ix],n[ix],false);
		for (it=0; it < (nt*ntr/n[ix]); it++) {
			i0 = sf_first_index (ix,it,1+1,n,s);
			for (iter=0; iter < 2; iter++) {
				sf_smooth2 (tr,i0,s[ix],false,window);
			}
		}
		sf_triangle_close(tr);
	}
	
	/* Tapering */
	pi = 4.0*atan(1.0);
	/*fprintf(stderr,"pi: %f\n",pi);
	fprintf(stderr,"ntr: %d\n",ntr);*/
	
	taper = (float *)calloc(ntr,sizeof(float));
	memset(taper,0,ntr*sizeof(float));

		
	for (ix=0; ix<151; ix++) {
		taper[ix] = (float)(0.5*(1.0-cos(2.0*pi*(ix-0.0)/300)));
		taper[ntr-ix-1] = taper[ix];
	}
	/*for (ix=(ntr-1); ix>(701-151-1); ix--) {
		taper[ix] = (float)(0.5*(1.0-cos(2.0*pi*(ix-0.0)/300)));

	}*/
	for (ix=151; ix<(ntr-151); ix++) {
		taper[ix] = 1.0;
	}
	/*for (ix=0; ix<ntr; ix++) {
		fprintf(stderr,"taper[%d]: %f\n",ix,taper[ix]);
	}*/
	
	FRefl = sf_input("refl");

	/*------------------------------------------------------------*/
	/* Loop over iterations */
	/*------------------------------------------------------------*/
	if (verb) fprintf(stderr,"Beginning of loop over iterations\n");
    for (iter=0; iter<niter; iter++) {
	
	/* Set Pminus and Qminus to 0 */
	memset(Pminus,0,2*nf*ntr*sizeof(float));
	memset(Qminus,0,2*nf*ntr*sizeof(float));

		/*------------------------------------------------------------*/
		/* Loop over shot positions */
		/*------------------------------------------------------------*/
		if (verb) fprintf(stderr,"Beginning of loop over shot positions\n");
		for (ishot=0; ishot<nshots; ishot++) {
   	
			/* Loop over receivers (traces) */
			#ifdef _OPENMP
			#pragma omp parallel for private(ix,it,a,b,c,d,e,f) \
				shared(Pminus,Qminus,Pplus,Qplus,taper,Refl)
			#endif 	
		  	for (ix=0; ix<ntr; ix++) {
				/* Loop over frequencies */
#if defined(__INTEL_COMPILER)
#pragma ivdep
#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC ivdep
#endif
				for (it=0; it<2*nf; it=it+2) {
					
					/*(x + yi)(u + vi) = (xu - yv) + (xv + yu)i*/
					a = Refl[ix*2*nf+it+ishot*2*nf*ntr]*taper[ishot];
					b = Refl[ix*2*nf+it+1+ishot*2*nf*ntr]*taper[ishot];
					c = Pplus[ishot*2*nf+it];
					d = Pplus[ishot*2*nf+it+1];
					e = Qplus[ishot*2*nf+it];
					f = Qplus[ishot*2*nf+it+1];

					Pminus[ix*2*nf+it]   += a*c - mode*b*d;
					Pminus[ix*2*nf+it+1] += mode*a*d + b*c;
					
					Qminus[ix*2*nf+it]   += a*e - mode*b*f;
					Qminus[ix*2*nf+it+1] += mode*a*f + b*e;
				
				} /* End of loop over frequencies */	
			} /* End of loop over receivers (traces) */
			
			if (verb) if(ishot%50==0) fprintf(stderr,"Trace %d\n",ishot);

		} /* End of loop over shot positions */

		/* Save a copy of pplus and qplus before creating their next iteration */
		memcpy(pplustemp,pplus,nt*ntr*sizeof(float));
		memcpy(qplustemp,qplus,nt*ntr*sizeof(float));

		/* Build the next iteration of Pplus and Qplus */
		fft1(Pminus,pminus,FRefl,1,0,1);
		fft1(Qminus,qminus,FRefl,1,0,1);
		
		if (verb) fprintf(stderr,"Build the next iteration of pplus and qplus\n");
		#ifdef _OPENMP
		#pragma omp parallel for private(ix,it) \
			shared(pminus,qminus,pplus,qplus,pplus0,window)
		#endif
		for (ix=0; ix<ntr; ix++) {
#if defined(__INTEL_COMPILER)
#pragma ivdep
#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC ivdep
#endif
			for (it=0; it<nt; it++) {
				pplus[it+ix*nt] = pplus0[it+ix*nt] - scale*window[it+ix*nt]*pminus[it+ix*nt];
				qplus[it+ix*nt] = pplus0[it+ix*nt] + scale*window[it+ix*nt]*qminus[it+ix*nt];
			}	
		}
		
		fft1(pplus,Pplus,Fplus,0,0,1);
		fft1(qplus,Qplus,Fplus,0,0,1);

		if (verb) fprintf(stderr,"%d %d\n",ix,it);

		if(iter%10==0) fprintf(stderr,"Iteration %d\n",iter);

	} /* End of loop over iterations */ 
	
	/* Build Gp and Gm */
	if (verb) fprintf(stderr,"Build Gp and Gm\n");
	#ifdef _OPENMP
	#pragma omp parallel for private(ix,it) \
		shared(Gp,Gm,G,H,pminus,qminus,pplustemp,qplustemp,pplus0)
	#endif
	for (ix=0; ix<ntr; ix++) {
#if defined(__INTEL_COMPILER)
#pragma ivdep
#elif defined(__GNUC__) && !defined(__clang__)
#pragma GCC ivdep
#endif
		for (it=0; it<nt; it++) {
			Gp[it+ix*nt] = 0.5*( pplustemp[it+ix*nt] + scale*pminus[it+ix*nt] + qplustemp[it+ix*nt] - scale*qminus[it+ix*nt] );
			Gm[it+ix*nt] = 0.5*( pplustemp[it+ix*nt] + scale*pminus[it+ix*nt] - qplustemp[it+ix*nt] + scale*qminus[it+ix*nt] );
			if (Gtot) {
				G[it+ix*nt] = pplustemp[it+ix*nt] + scale*pminus[it+ix*nt];
			}
			if (Htot) {
				H[it+ix*nt] = qplustemp[it+ix*nt] - scale*qminus[it+ix*nt];
			}
			/*p[it+ix*nt] = 1.0*( pplus[it+ix*nt] + scale*pminus[it+ix*nt] );
			q[it+ix*nt] = 1.0*( qplus[it+ix*nt] - scale*qminus[it+ix*nt] );*/
		}	
	}

	/* Write the final result */
    /*FRefl = sf_input(argv[1]);*/
	/*fft1(Gp,,FRefl,1,0,0);
	fft1(Gm,,FRefl,1,0,0);*/
	
	sf_floatwrite(Gp,nt*ntr,FGp);
	sf_floatwrite(Gm,nt*ntr,FGm);
	if (Gtot) {
		sf_floatwrite(G,nt*ntr,FG);
		sf_fileclose(FG);
		free(G);
	}
	
	if (Htot) {
		sf_floatwrite(H,nt*ntr,FH);
		sf_fileclose(FH);
		free(H);
	}

	if (pandq) {
		sf_floatwrite(pplustemp,nt*ntr,Fp);
		sf_floatwrite(pminus,nt*ntr,Fq);
		sf_fileclose(Fp);
		sf_fileclose(Fq);
	}
	
	if (twin) {
		sf_floatwrite(window,nt*ntr,Ftwin);
		sf_fileclose(Ftwin);
	}
	sf_fileclose(Fplus);
	sf_fileclose(FRefl);
	sf_fileclose(FGp);
	sf_fileclose(FGm);
	
	free(Gp);
	free(Gm);
	free(pplus);
	free(pplusinv);
	free(pplustemp);
	free(Pplus);
	free(Pplus_trace);
	free(Pminus);
	free(qplus);
	free(qplustemp);
	free(Qplus);
	free(Qplus_trace);
	free(Qminus);
	free(Refl);
	free(window);
	free(tw);
	free(filename1);
	
    exit (0);
}

