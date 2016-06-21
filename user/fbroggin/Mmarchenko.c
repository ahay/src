/* Marchenko-Wapenaar-Broggini iterative scheme

sfmarchenko < downgoing.rsf refl=REFL_000.rsf conj=y verb=n Gtot=y niter=21 nshots=401 scale=1 eps=1e-4 shift=5 Gm=Gm.rsf G=G.rsf> Gp.rsf

======= INPUTS ============

p0plus.rsf: initial downgoing wave field

REFL_000.rsf: Fourier transform of the reflection response

======= PARAMETERS ========

conj  = [y]/n	- complex-conjugation of the first input (corresponds to time-reversal in time)
verb = y/[n]	- verbosity flag
twin  = y/[n]	- returns the timewindow as one of the outputs (window=window.rsf)
pandq  = y/[n]	- pandq=true: returns p and q, pandq=false returns Gp and Gm
Gtot  = y/[n]	- Gtot=true returns G=Gp+Gm
Htot  = y/[n]	- Htot=true returns H=Gp-Gm
niter  = 1		- number of iterations
ntaper  = 101	- tapering width for each side
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

void fft1(float *, float *, sf_file, bool, bool, bool);

int main(int argc, char* argv[]) {

	bool verb, conj, twin, pandq, Gtot, Htot;

	/* OMP parameters */
#ifdef _OPENMP
	int ompnth;
#endif 	

	float *pplus0, *pplus, *pplusinv, *pplustemp, *Pplus, *Pplus_trace, *pminus, *Pminus, *Refl, *Gp, *Gm, *G, *H;
	float *qplus, *qplusinv, *qplustemp, *Qplus, *Qplus_trace, *qminus, *Qminus;
	float *window, *window_all, *taper, pi;
	int *tw;

	/* I/O files */
	sf_file Fplus;
	sf_file FRefl;
	sf_file FGp;
	sf_file FGm;
	sf_file FG;
	sf_file FH;
	sf_file Ftwin;
	sf_file Fp;
	sf_file Fq;

	/* Cube axes */
	sf_axis at, af, ax, as, avs;

	int nt, nf, ntr, nshots, nvs, mode, niter, len, ntaper;
	int i, it, ix, ishot, iter, ivs, i0;
	int twc, twa, shift, n[2], rect[2], s[2];
	float scale, eps, dt, df, dx, ds, dvs, ot, of, a, b, c, d, e, f;

	sf_triangle tr;

	/*------------------------------------------------------------*/
	/* Initialize RSF parameters 								  */
	/*------------------------------------------------------------*/
	sf_init(argc, argv);

	/*------------------------------------------------------------*/
	/* Initialize OMP parameters */
	/*------------------------------------------------------------*/
#ifdef _OPENMP
	ompnth = omp_init();
#endif	

	/*------------------------------------------------------------*/
	/* Flags 													  */
	/*------------------------------------------------------------*/
	if (!sf_getbool("verb", &verb)) verb = false; /* verbosity flag */
	if (!sf_getbool("conj", &conj)) conj = false; /* complex conjugation (time-reversal) flag */
	if (!sf_getbool("twin", &twin)) twin = false; /* returns the timewindow as one of the outputs */
	if (!sf_getbool("pandq", &pandq)) pandq = false; /* pandq=true: returns p and q */
	if (!sf_getbool("Gtot", &Gtot)) Gtot = false; /* Gtot=true: returns G=Gp+Gm */
	if (!sf_getbool("Htot", &Htot)) Htot = false; /* Htot=true: returns H=Gp-Gm */
	if (!sf_getint("niter", &niter)) niter = 1; /* number of iterations */
	if (!sf_getint("ntaper", &ntaper)) ntaper = 101; /* tapering widht */
	if (!sf_getfloat("scale", &scale)) scale = 1.0; /* scale factor */
	if (!sf_getfloat("eps", &eps)) eps = 1e-4; /* threshold for the timewindow */
	if (!sf_getint("shift", &shift)) shift = 5; /* shift in samples for the timewindow */

	if (verb) {
		fprintf(stderr, "This program was called with \"%s\".\n", argv[0]);
		/*fprintf(stderr,"Nr: %d Nx: %d Nt:%d\n",nr,nx,nt);*/

		if (argc > 1) {
			for (i = 1; i < argc; i++) {
				fprintf(stderr, "argv[%d] = %s\n", i, argv[i]);
			}
		} else {
			fprintf(stderr, "The command had no other arguments.\n");
		}
	}
	
	/* Check if ntaper is odd */
	if (!(ntaper % 2)) {
		ntaper += 1;
		fprintf(stderr, "ntaper is %d\n", ntaper);
	}

	/*------------------------------------------------------------*/
	/* I/O files 												  */
	/*------------------------------------------------------------*/
	/* "in" is the transposed version of p00plus_xxxx.rsf
	   where xxxx denotes the depth of the virtual sources
	   Dimensions of p00plus_xxxx.rsf should be: n1=nt,n2=ntr,n3=nvs */
	Fplus = sf_input("in");

	/* refl is refl_fft_all.rsf
	   It is used to read nf, df, of
	   Dimensions are: n1=nf,n2=ntr,n3=nshots */
	FRefl = sf_input("refl");

	FGp = sf_output("out"); /* Gp */
	FGm = sf_output("Gm");

	if (Gtot) {
		FG = sf_output("G");
	}

	if (Htot) {
		FH = sf_output("H");
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
	/* Time */
	at = sf_iaxa(Fplus, 1);	sf_setlabel(at, "Time"); if (verb) sf_raxa(at);
	/* Space - receivers at the acquisition surface */
	ax = sf_iaxa(Fplus, 2); sf_setlabel(ax, "r"); if (verb) sf_raxa(ax);
	/* Virtual sources */
	avs = sf_iaxa(Fplus, 3); sf_setlabel(avs, "VS"); if (verb) sf_raxa(avs);
	
	/* Frequency */
	af = sf_iaxa(FRefl, 1); sf_setlabel(af, "Frequency"); if (verb) sf_raxa(af);
	/* Shots */
	as = sf_iaxa(FRefl, 3); sf_setlabel(as, "Shots"); if (verb) sf_raxa(as);
	
	/* Number of time steps and time sampling */
	nt = sf_n(at);
	dt = sf_d(at);
	ot = sf_o(at);
	/* Number of frequencies and frequency sampling */
	nf = sf_n(af);
	df = sf_d(af);
	of = sf_o(af);
	/* Number of traces and spacing */
	ntr = sf_n(ax);
	dx = sf_d(ax);
	/* Number of shots and spacing */
	nshots = sf_n(as);
	ds = sf_d(as);
	/* Number of virtual sources and spacing */
	nvs = sf_n(avs);
	dvs = sf_d(avs);

	if (verb) fprintf(stderr, "nt: %d nf: %d ntr: %d nshots: %d nvs: %d\n", nt, nf, ntr, nshots, nvs);

	/*------------------------------------------------------------*/
	/* Setup output data and wave field header					  */
	/*------------------------------------------------------------*/
	sf_oaxa(FGp, at, 1);
	sf_oaxa(FGp, ax, 2);
	sf_oaxa(FGp, avs, 3);
	
	sf_oaxa(FGm, at, 1);
	sf_oaxa(FGm, ax, 2);
	sf_oaxa(FGm, avs, 3);
	
	if (Gtot) {
		sf_oaxa(FG, at, 1);
		sf_oaxa(FG, ax, 2);
		sf_oaxa(FG, avs, 3);
	}

	if (Htot) {
		sf_oaxa(FH, at, 1);
		sf_oaxa(FH, ax, 2);
		sf_oaxa(FH, avs, 3);
	}

	if (pandq) {
		sf_oaxa(Fp, at, 1);
		sf_oaxa(Fp, ax, 2);
		sf_oaxa(Fp, avs, 3);
		sf_oaxa(Fq, at, 1);
		sf_oaxa(Fq, ax, 2);
		sf_oaxa(Fq, avs, 3);
	}

	if (twin) {
		sf_oaxa(Ftwin, at, 1);
		sf_oaxa(Ftwin, ax, 2);
		sf_oaxa(Ftwin, avs, 3);
	}

	/*------------------------------------------------------------*/
	/* Allocate arrays											  */
	/*------------------------------------------------------------*/
	
	/* Downgoing wave fields - Time */
	pplus0 = (float *) calloc(nt*ntr, sizeof (float));
	pplus = (float *) calloc(nt*ntr, sizeof (float));
	pplustemp = (float *) calloc(nt*ntr, sizeof (float));
	/*pplusinv = (float *) calloc(nt*ntr, sizeof (float));*/
	qplus = (float *) calloc(nt*ntr, sizeof (float));
	qplustemp = (float *) calloc(nt*ntr, sizeof (float));

	/* Downgoing wave fields - Frequency */
	Pplus = (float *) calloc(2 * nf*ntr, sizeof (float));
	Qplus = (float *) calloc(2 * nf*ntr, sizeof (float));

	/* Upgoing wave fields - Time */
	pminus = (float *) calloc(nt*ntr, sizeof (float));
	qminus = (float *) calloc(nt*ntr, sizeof (float));

	/* Downgoing wave fields - Frequency */
	Pminus = (float *) calloc(2 * nf*ntr, sizeof (float));
	Qminus = (float *) calloc(2 * nf*ntr, sizeof (float));

	/* This is currently NOT used */
	/* Transpose of p00plus_xxxx_xxxx */
	/*for (ix = 0; ix < ntr; ix++) {
		for (it = 0; it < nt; it++) {
			pplusinv[ix * ntr + it] = pplus0[it * ntr + ix];
		}
	}*/

	/* Single trace (in frequency) of the downgoing wave field */
	Pplus_trace = (float *) calloc(2 * nf, sizeof (float));
	Qplus_trace = (float *) calloc(2 * nf, sizeof (float));
	
	/* Time window */
	tw = (int *) calloc(ntr, sizeof (int));
	window = (float *) calloc(nt*ntr, sizeof (float));
	if (twin) {
		window_all = (float *) calloc(nt*ntr*nvs, sizeof (float));
	}

	/* Output wave fields */
	Gp = (float *) calloc(nt*ntr*nvs, sizeof (float));
	Gm = (float *) calloc(nt*ntr*nvs, sizeof (float));
	if (Gtot) {
		G = (float *) calloc(nt*ntr*nvs, sizeof (float));
	}

	if (Htot) {
		H = (float *) calloc(nt*ntr*nvs, sizeof (float));
	}

	/* Time-reversal flag */
	if (conj) {
		mode = -1;
	} else {
		mode = +1;
	}

	/* Load the reflection response into the memory */
	if (verb) fprintf(stderr, "Before loading R %d\n", 2 * nf * ntr);
	Refl = (float *) calloc(2 * nf * ntr * nshots, sizeof (float));
	sf_floatread(&Refl[0], 2 * nf * ntr * nshots, FRefl);

	/* Tapering */
	pi = 4.0 * atan(1.0);
	
	taper = (float *) calloc(ntr, sizeof (float));
	memset(taper, 0, ntr * sizeof (float));

	/* Left and right tapering */
	for (ix = 0; ix < ntaper; ix++) {
		taper[ix] = (float) (0.5 * (1.0 - cos(2.0 * pi * (ix - 0.0) / ((ntaper - 1)*2))));
		taper[ntr - ix - 1] = taper[ix];
	}
	/* Central part of the taper is equal to 1 */
	for (ix = ntaper; ix < (ntr - ntaper); ix++) {
		taper[ix] = 1.0;
	}
	/* Print taper values */
	/*for (ix = 0; ix < ntr; ix++) {
		fprintf(stderr, "taper[%d]: %f\n", ix, taper[ix]);
	}*/

	/*------------------------------------------------------------*/
	/* Loop over VS */
	/*------------------------------------------------------------*/
	if (verb) fprintf(stderr, "Beginning of loop over VS\n");
	for (ivs = 0; ivs < nvs; ivs++) {
		
		if (ivs % 1 == 0) fprintf(stderr, "VS %d\n", ivs);

		/* Read first arrival */
		sf_seek(Fplus, ivs * nt * ntr * sizeof (float), SEEK_SET);
		sf_floatread(pplus0, nt * ntr, Fplus);
		memcpy(pplus, pplus0, nt * ntr * sizeof (float));

		/* The three flags of fft1 are: inv, sym, and opt */
		fft1(pplus0, Pplus, Fplus, 0, 0, 1);
		memcpy(Qplus, Pplus, 2 * nf * ntr * sizeof (float));

		/* Build time-window */
		memset(window, 0, nt * ntr * sizeof (float));
		
		if (verb) fprintf(stderr, "Build the time-window\n");
		#ifdef _OPENMP
		#pragma omp parallel for private(ix,it) \
			shared(pplus0,tw)
		#endif
		for (ix = 0; ix < ntr; ix++) {
			for (it = 0; it < nt; it++) {
				if (pplus0[it + ix * nt] > eps) {
					/*tw[ix] = it*dt+ot;*/
					tw[ix] = it;
					/*if (verb) fprintf(stderr,"%d %d\n",ix,it);*/
					break;
				}
			}
		}
		#ifdef _OPENMP
		#pragma omp parallel for private(ix,it,twa,twc) \
			shared(window)
		#endif
		for (ix = 0; ix < ntr; ix++) {
			twc = (int) (tw[ix] - shift);
			twa = (int) (-twc + shift + nt);
			/*if (verb) fprintf(stderr,"%d %d\n",twc,twa);*/
			for (it = 0; it < nt; it++) {
				if ((it < twc) || (it > twa)) {
					window[it + ix * nt] = 1.0;
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

		for (ix = 0; ix <= 1; ix++) {
			if (rect[ix] <= 1) continue;
			tr = sf_triangle_init(rect[ix], n[ix],false);
			for (it = 0; it < (nt * ntr / n[ix]); it++) {
				i0 = sf_first_index(ix, it, 1 + 1, n, s);
				for (iter = 0; iter < 2; iter++) {
					sf_smooth2(tr, i0, s[ix], false, window);
				}
			}
			sf_triangle_close(tr);
		}
		
		if (twin) {
			memcpy(&window_all[ivs * nt * ntr], window, nt * ntr * sizeof (float));
		}

		/*------------------------------------------------------------*/
		/* Loop over iterations */
		/*------------------------------------------------------------*/
		if (verb) fprintf(stderr, "Beginning of loop over iterations\n");
		for (iter = 0; iter < niter; iter++) {

			/* Set Pminus and Qminus to 0 */
			memset(Pminus, 0, 2 * nf * ntr * sizeof (float));
			memset(Qminus, 0, 2 * nf * ntr * sizeof (float));

			/*------------------------------------------------------------*/
			/* Loop over shot positions */
			/*------------------------------------------------------------*/
			if (verb) fprintf(stderr, "Beginning of loop over shot positions\n");
			for (ishot = 0; ishot < nshots; ishot++) {

				/*------------------------------------------------------------*/
				/* Loop over receivers (traces) */
				/*------------------------------------------------------------*/
				#ifdef _OPENMP
				#pragma omp parallel for private(ix,it,a,b,c,d,e,f) \
					shared(Pminus,Qminus,Pplus,Qplus,taper,Refl,ivs)
				#endif 	
				for (ix = 0; ix < ntr; ix++) {
					/*------------------------------------------------------------*/
					/* Loop over frequencies */
					/*------------------------------------------------------------*/
					#pragma ivdep
					for (it = 0; it < 2 * nf; it = it + 2) {

						/*(x + yi)(u + vi) = (xu - yv) + (xv + yu)i*/
						a = Refl[ix * 2 * nf + it + ishot * 2 * nf * ntr] * taper[ishot];
						b = Refl[ix * 2 * nf + it + 1 + ishot * 2 * nf * ntr] * taper[ishot];
						c = Pplus[ishot * 2 * nf + it];
						d = Pplus[ishot * 2 * nf + it + 1];
						e = Qplus[ishot * 2 * nf + it];
						f = Qplus[ishot * 2 * nf + it + 1];

						Pminus[ix * 2 * nf + it] += a * c - mode * b*d;
						Pminus[ix * 2 * nf + it + 1] += mode * a * d + b*c;

						Qminus[ix * 2 * nf + it] += a * e - mode * b*f;
						Qminus[ix * 2 * nf + it + 1] += mode * a * f + b*e;

					} /* End of loop over frequencies */
				} /* End of loop over receivers (traces) */

				if (verb) if (ishot % 50 == 0) fprintf(stderr, "Trace %d\n", ishot);

			} /* End of loop over shot positions */

			/* Save a copy of pplus and qplus before creating their next iteration */
			memcpy(pplustemp, pplus, nt * ntr * sizeof (float));
			memcpy(qplustemp, qplus, nt * ntr * sizeof (float));

			/* Build the next iteration of Pplus and Qplus */
			fft1(Pminus, pminus, FRefl, 1, 0, 1);
			fft1(Qminus, qminus, FRefl, 1, 0, 1);

			if (verb) fprintf(stderr, "Build the next iteration of pplus and qplus\n");
			#ifdef _OPENMP
			#pragma omp parallel for private(ix,it) \
				shared(pminus,qminus,pplus,qplus,pplus0,window)
			#endif
			for (ix = 0; ix < ntr; ix++) {
				#pragma ivdep
				for (it = 0; it < nt; it++) {
					pplus[it + ix * nt] = pplus0[it + ix * nt] - scale * window[it + ix * nt] * pminus[it + ix * nt];
					qplus[it + ix * nt] = pplus0[it + ix * nt] + scale * window[it + ix * nt] * qminus[it + ix * nt];
				}
			}

			fft1(pplus, Pplus, Fplus, 0, 0, 1);
			fft1(qplus, Qplus, Fplus, 0, 0, 1);

			if (verb) fprintf(stderr, "%d %d\n", ix, it);

			if (iter % 10 == 0) fprintf(stderr, "Iteration %d\n", iter);
			
		} /* End of loop over iterations */

		/* Build Gp and Gm */
		if (verb) fprintf(stderr, "Build Gp and Gm\n");
		#ifdef _OPENMP
		#pragma omp parallel for private(ix,it) \
			shared(Gp,Gm,G,H,pminus,qminus,pplustemp,qplustemp,pplus0,ivs)
		#endif
		for (ix = 0; ix < ntr; ix++) {
			#pragma ivdep
			for (it = 0; it < nt; it++) {
				Gp[it + ix * nt + ivs * ntr * nt] = 0.5 * (pplustemp[it + ix * nt] + scale * pminus[it + ix * nt] + qplustemp[it + ix * nt] - scale * qminus[it + ix * nt]);
				Gm[it + ix * nt + ivs * ntr * nt] = 0.5 * (pplustemp[it + ix * nt] + scale * pminus[it + ix * nt] - qplustemp[it + ix * nt] + scale * qminus[it + ix * nt]);
				if (Gtot) {
					G[it + ix * nt + ivs * ntr * nt] = pplustemp[it + ix * nt] + scale * pminus[it + ix * nt];
				}
				if (Htot) {
					H[it + ix * nt + ivs * ntr * nt] = qplustemp[it + ix * nt] - scale * qminus[it + ix * nt];
				}
				/*p[it+ix*nt] = 1.0*( pplus[it+ix*nt] + scale*pminus[it+ix*nt] );
				q[it+ix*nt] = 1.0*( qplus[it+ix*nt] - scale*qminus[it+ix*nt] );*/
			}
		}
	} /* End of loop over VS */

	/* Write the final result */
	/*FRefl = sf_input(argv[1]);*/
	/*fft1(Gp,,FRefl,1,0,0);
	fft1(Gm,,FRefl,1,0,0);*/
	
	fprintf(stderr, "Break 1\n");
	sf_floatwrite(Gp, nt*ntr*nvs, FGp);
	fprintf(stderr, "Break 2\n");
	sf_floatwrite(Gm, nt*ntr*nvs, FGm);
	fprintf(stderr, "Break 3\n");
	if (Gtot) {
		sf_floatwrite(G, nt*ntr*nvs, FG);
		sf_fileclose(FG);
		free(G);
	}
	fprintf(stderr, "Break 4\n");

	if (Htot) {
		sf_floatwrite(H, nt*ntr*nvs, FH);
		sf_fileclose(FH);
		free(H);
	}
	fprintf(stderr, "Break 5\n");

	if (pandq) {
		sf_floatwrite(pplustemp, nt*ntr, Fp);
		sf_floatwrite(pminus, nt*ntr, Fq);
		sf_fileclose(Fp);
		sf_fileclose(Fq);
	}
	fprintf(stderr, "Break 6\n");

	if (twin) {
		sf_floatwrite(window_all, nt*ntr*nvs, Ftwin);
		sf_fileclose(Ftwin);
		free(window_all);
	}
	fprintf(stderr, "Break 7\n");
	sf_fileclose(Fplus);
	sf_fileclose(FRefl);
	sf_fileclose(FGp);
	sf_fileclose(FGm);
	fprintf(stderr, "Break 8\n");

	free(Gp);
	free(Gm);
	free(pplus);
	/*free(pplusinv);*/
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
	fprintf(stderr, "Break 9\n");

	exit(0);
}

