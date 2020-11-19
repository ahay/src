/* 2D Helmholtz solver by LU factorization. */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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
#include <umfpack.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fdprep.h"

int main(int argc, char* argv[])
{
    bool verb, save, load, hermite;
    int n1, n2, npml, pad1, pad2, is, ns, iw, nw;
    SuiteSparse_long n, nz, *Ti, *Tj;
    float d1, d2, **v, ds, os, dw, ow;
    double omega, *Tx, *Tz;
    SuiteSparse_long *Ap, *Ai, *Map;
    double *Ax, *Az, **Xx, **Xz, **Bx, **Bz;
    void *Symbolic, **Numeric;
    double Control[UMFPACK_CONTROL];
    sf_complex ***f;
    char *datapath, *insert, *append;
    size_t srclen, inslen;
    sf_file in, out, source;
    int uts, its, mts;
    sf_timer timer;
    char *order;

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */
    
    if (verb)
	timer = sf_timer_init();
    else
	timer = NULL;

    if (!sf_getbool("save",&save)) save=false;
    /* save LU */

    if (!sf_getbool("load",&load)) load=false;
    /* load LU */

    if (save || load) {
	datapath = sf_histstring(in,"in");
	srclen = strlen(datapath);
	insert = sf_charalloc(6);
    } else {
	datapath = NULL;
	srclen = 0;
	insert = NULL;
	append = NULL;
    }

    if (!sf_getint("uts",&uts)) uts=0;
    /* number of OMP threads */

#ifdef _OPENMP
    mts = omp_get_max_threads();
#else
    mts = 1;
#endif

    uts = (uts < 1)? mts: uts;
    if (verb) sf_warning("Using %d out of %d threads.",uts,mts);

    if (!sf_getbool("hermite",&hermite)) hermite=false;
    /* Hermite operator */
    
    if (!sf_getint("npml",&npml)) npml=10;
    /* PML width */

    if (NULL == (order = sf_getstring("order"))) order="j";
    /* discretization scheme (default optimal 9-point) */

    fdprep_order(order);

    /* read input dimension */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input.");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input.");

    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input.");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input.");

    /* read input */
    v = sf_floatalloc2(n1,n2);
    sf_floatread(v[0],n1*n2,in);
    
    /* read source */
    if (NULL == sf_getstring("source"))
	sf_error("Need source=");
    source = sf_input("source");

    if (!sf_histint(source,"n3",&ns)) sf_error("No ns=.");
    if (!sf_histfloat(source,"d3",&ds)) ds=d2;
    if (!sf_histfloat(source,"o3",&os)) os=0.;
    if (!sf_histint(source,"n4",&nw)) sf_error("No nw=.");
    if (!sf_histfloat(source,"d4",&dw)) sf_error("No dw=.");
    if (!sf_histfloat(source,"o4",&ow)) sf_error("No ow=.");

    f = sf_complexalloc3(n1,n2,ns);

    /* write output header */
    sf_settype(out,SF_COMPLEX);

    sf_putint(out,"n3",ns);
    sf_putfloat(out,"d3",ds);
    sf_putfloat(out,"o3",os);
    sf_putstring(out,"label3","Shot");
    sf_putstring(out,"unit3","");
    sf_putint(out,"n4",nw);
    sf_putfloat(out,"d4",dw);
    sf_putfloat(out,"o4",ow);
    sf_putstring(out,"label4","Frequency");
    sf_putstring(out,"unit4","Hz");

    /* allocate temporary memory */
    if (load) {
	Ti = NULL; Tj = NULL; Tx = NULL; Tz = NULL; 
	Ap = NULL; Ai = NULL; Map = NULL; Ax = NULL; Az = NULL;
    }

    Bx = (double**) sf_alloc(uts,sizeof(double*));
    Bz = (double**) sf_alloc(uts,sizeof(double*));
    Xx = (double**) sf_alloc(uts,sizeof(double*));
    Xz = (double**) sf_alloc(uts,sizeof(double*));

    Numeric = (void**) sf_alloc(uts,sizeof(void*));

    /* LU control */
    umfpack_zl_defaults (Control);
    Control [UMFPACK_IRSTEP] = 0;

    /* loop over frequency */
    for (iw=0; iw < nw; iw++) {
	omega = (double) 2.*SF_PI*(ow+iw*dw);

	/* PML padding */
	pad1 = n1+2*npml;
	pad2 = n2+2*npml;

	n  = fdprep_n (pad1,pad2);
	nz = fdprep_nz(pad1,pad2);

	if (!load) {
	    Ti = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long));
	    Tj = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long));
	    Tx = (double*) sf_alloc(nz,sizeof(double));
	    Tz = (double*) sf_alloc(nz,sizeof(double));
	    
	    Ap = (SuiteSparse_long*) sf_alloc(n+1,sizeof(SuiteSparse_long));
	    Ai = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long));
	    Map = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long));
	    
	    Ax = (double*) sf_alloc(nz,sizeof(double));
	    Az = (double*) sf_alloc(nz,sizeof(double));
	}
	for (its=0; its < uts; its++) {
	    Bx[its] = (double*) sf_alloc(n,sizeof(double));
	    Bz[its] = (double*) sf_alloc(n,sizeof(double));
	    Xx[its] = (double*) sf_alloc(n,sizeof(double));
	    Xz[its] = (double*) sf_alloc(n,sizeof(double));
	}

	if (verb) {
	    sf_warning("Frequency %d of %d. Pad %d for PML.",iw+1,nw,npml);
	    sf_timer_start(timer);
	}

	/* LU file (append _lu* after velocity file) */
	if (save || load) {
	    sprintf(insert,"_lu%d",iw);
	    inslen = strlen(insert);
	    
	    append = malloc(srclen+inslen+1);

	    memcpy(append,datapath,srclen-5);
	    memcpy(append+srclen-5,insert,inslen);
	    memcpy(append+srclen-5+inslen,datapath+srclen-5,5+1);
	}

	if (!load) {
	    /* assemble matrix */	    
	    fdprep(omega,
		   n1, n2, d1, d2, v,
		   npml, pad1, pad2,
		   Ti, Tj, Tx, Tz);
	    
	    umfpack_zl_triplet_to_col (n, n, nz, 
			Ti, Tj, Tx, Tz, 
			Ap, Ai, Ax, Az, Map);
	    
	    /* LU */
	    umfpack_zl_symbolic (n, n, 
			Ap, Ai, Ax, Az, 
			&Symbolic, Control, NULL);

	    umfpack_zl_numeric (Ap, Ai, Ax, Az, 
	        Symbolic, &Numeric[0], 
		    Control, NULL);

	    /* save Numeric */
#ifdef _OPENMP
	    umfpack_zl_save_numeric (Numeric[0], append);
	    
	    for (its=1; its < uts; its++) {
		umfpack_zl_load_numeric (&Numeric[its], append);
	    }
	    
	    if (!save) {
		(void) remove (append);
		(void) remove ("numeric.umf");
	    }
#else
	    if (save) umfpack_zl_save_numeric (Numeric[0], append);
#endif
	} else {
	    /* load Numeric */
	    for (its=0; its < uts; its++) {
		umfpack_zl_load_numeric (&Numeric[its], append);
	    }
	}

	if (save || load) free(append);

	/* read source */
	sf_complexread(f[0][0],n1*n2*ns,source);

	/* loop over shots */
#ifdef _OPENMP
#pragma omp parallel for num_threads(uts) private(its)
#endif
	for (is=0; is < ns; is++) {
#ifdef _OPENMP
	    its = omp_get_thread_num();
#else
	    its = 0;
#endif
	    
	    fdpad(npml,pad1,pad2, f[is],Bx[its],Bz[its]);

	    umfpack_zl_solve (hermite? UMFPACK_At: UMFPACK_A, 
		    NULL, NULL, NULL, NULL, 
			Xx[its], Xz[its], Bx[its], Bz[its], 
			Numeric[its], Control, NULL);	    
	    	    
	    fdcut(npml,pad1,pad2, f[is],Xx[its],Xz[its]);
	}

	/* write wavefields */
	sf_complexwrite(f[0][0],n1*n2*ns,out);

	if (verb) {
	    sf_timer_stop (timer);
	    sf_warning("Finished in %g seconds.",sf_timer_get_diff_time(timer)/1.e3);
	}

	if (!load) (void) umfpack_zl_free_symbolic (&Symbolic);
	for (its=0; its < uts; its++) {
	    (void) umfpack_zl_free_numeric (&Numeric[its]);
	}

	if (!load) {
	    free(Ti); free(Tj); free(Tx); free(Tz);
	    free(Ap); free(Ai); free(Map);
	    free(Ax); free(Az);
	}
	for (its=0; its < uts; its++) {
	    free(Bx[its]); free(Bz[its]); free(Xx[its]); free(Xz[its]);
	}
    }

    exit(0);
}
