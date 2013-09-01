/* IWI interface for ODCIG migration */
/*
  Copyright (C) 2013 University of Texas at Austin
  
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
#include "iwimodl.h"

static float d1, d2, ow, dw, ****image;
static int n1, n2, nh, ns, nw, npml;
static sf_file source, data;
static sf_fslice sfile, rfile;
static SuiteSparse_long *Ti, *Tj;
static SuiteSparse_long *Ap, *Ai, *Map;
static void *Symbolic, **Numeric;
static double Control[UMFPACK_CONTROL];
static double *Tx, *Tz;
static double *Ax, *Az, **Xx, **Xz, **Bx, **Bz;
static sf_complex ***srce, ***recv;
static bool load;
static int uts;
static char *datapath, *insert, *append;
static size_t srclen, inslen;

void iwimodl_init(char *order,
		  int npml0,
		  int nn1, int nn2, 
		  float dd1, float dd2,
		  int nh0, int ns0, 
		  float ow0, float dw0, int nw0,
		  sf_file source0, sf_file data0,
		  sf_fslice sfile0, sf_fslice rfile0,
		  bool load0, char *datapath0,
		  int uts0)
/*< initialization >*/
{
    int its;

    fdprep_order(order);

    npml = npml0;
    
    n1 = nn1;
    n2 = nn2;
    d1 = dd1;
    d2 = dd2;

    nh = nh0; ns = ns0;
    ow = ow0; dw = dw0; nw = nw0;

    source = source0;
    data = data0;

    sfile = sfile0;
    rfile = rfile0;

    load = load0;
    datapath = datapath0;
    
    uts = uts0;

    /* LU file */
    if (load) {
	srclen = strlen(datapath);
	insert = sf_charalloc(6);
    } else {
	srclen = 0;
	insert = NULL;
	append = NULL;
    }    
    
    /* allocate temporary space */
    srce = sf_complexalloc3(n1,n2,ns);
    recv = sf_complexalloc3(n1,n2,ns);
    
    Bx = (double**) sf_alloc(uts,sizeof(double*));
    Bz = (double**) sf_alloc(uts,sizeof(double*));
    Xx = (double**) sf_alloc(uts,sizeof(double*));
    Xz = (double**) sf_alloc(uts,sizeof(double*));

    image = (float****) sf_alloc(uts,sizeof(float***));
    for (its=0; its < uts; its++) {
	image[its] = sf_floatalloc3(n1,n2,2*nh+1);
    }

    Numeric = (void**) sf_alloc(uts,sizeof(void*));

    /* LU control */
    umfpack_zl_defaults (Control);
    Control [UMFPACK_IRSTEP] = 0;
}

void iwimodl_free()
/*< free >*/
{
    int its;

    free(srce[0][0]); free(srce[0]); free(srce);
    free(recv[0][0]); free(recv[0]); free(recv);

    for (its=0; its < uts; its++) {
	free(image[its][0][0]); free(image[its][0]); free(image[its]);
    }
    free(image);
}

void iwimodl_clean()
/*< clean up >*/
{
    int its, ih, j, i;

#ifdef _OPENMP
#pragma omp parallel for num_threads(uts) private(ih,j,i)
#endif
    for (its=0; its < uts; its++) {
	for (ih=-nh; ih < nh+1; ih++) {
	    for (j=0; j < n2; j++) {
		for (i=0; i < n1; i++) {
		    image[its][ih+nh][j][i] = 0.;
		}
	    }
	}
    }
}

void iwimodl_modl(float **vel   /* current velocity */,
		  float *image0 /* extended image */)
/*< modeling >*/
{
    int iw, is, its, ih, j, i;
    int pad1, pad2;
    SuiteSparse_long n, nz;
    double omega;

    /* PML padding */
    pad1 = n1+2*npml;
    pad2 = n2+2*npml;

    n  = fdprep_n (pad1,pad2);
    nz = fdprep_nz(pad1,pad2);

    /* loop over frequency */
    for (iw=0; iw < nw; iw++) {
	omega = (double) 2.*SF_PI*(ow+iw*dw);	

	Ti = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long));
	Tj = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long));
	Tx = (double*) sf_alloc(nz,sizeof(double));
	Tz = (double*) sf_alloc(nz,sizeof(double));
	    
	Ap = (SuiteSparse_long*) sf_alloc(n+1,sizeof(SuiteSparse_long));
	Ai = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long));
	Map = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long));

	Ax = (double*) sf_alloc(nz,sizeof(double));
	Az = (double*) sf_alloc(nz,sizeof(double));

	for (its=0; its < uts; its++) {
	    Bx[its] = (double*) sf_alloc(n,sizeof(double));
	    Bz[its] = (double*) sf_alloc(n,sizeof(double));
	    Xx[its] = (double*) sf_alloc(n,sizeof(double));
	    Xz[its] = (double*) sf_alloc(n,sizeof(double));
	}

	/* LU file (append _inv* after velocity file) */
	if (load) {
	    sprintf(insert,"_inv%d",iw);
	    inslen = strlen(insert);
	    
	    append = malloc(srclen+inslen+1);

	    memcpy(append,datapath,srclen-5);
	    memcpy(append+srclen-5,insert,inslen);
	    memcpy(append+srclen-5+inslen,datapath+srclen-5,5+1);
	}

	/* assemble matrix */
	fdprep(omega,
	       n1, n2, d1, d2, vel,
	       npml, pad1, pad2,
	       Ti, Tj, Tx, Tz);	    
	
	(void) umfpack_zl_triplet_to_col (n, n, nz, 
					  Ti, Tj, Tx, Tz, 
					  Ap, Ai, Ax, Az, Map);	    

	/* LU */
	(void) umfpack_zl_symbolic (n, n, 
				    Ap, Ai, Ax, Az, 
				    &Symbolic, Control, NULL);
	
	(void) umfpack_zl_numeric (Ap, Ai, Ax, Az, 
				   Symbolic, &Numeric[0], 
				   Control, NULL);
	
	/* save Numeric */
#ifdef _OPENMP
	(void) umfpack_zl_save_numeric (Numeric[0], append);
	
	for (its=1; its < uts; its++) {
	    (void) umfpack_zl_load_numeric (&Numeric[its], append);
	}
	
	if (!load) {
	    (void) remove (append);
	    (void) remove ("numeric.umf");
	}
#else
	if (load) (void) umfpack_zl_save_numeric (Numeric[0], append);
#endif	
	
	if (load) free(append);

	/* read source and data */
	sf_complexread(srce[0][0],n1*n2*ns,source);
	sf_complexread(recv[0][0],n1*n2*ns,data);

	/* loop over shots */
#ifdef _OPENMP
#pragma omp parallel for num_threads(uts) private(its,ih,j,i)
#endif
	for (is=0; is < ns; is++) {
#ifdef _OPENMP
	    its = omp_get_thread_num();
#else
	    its = 0;
#endif

	    /* source wavefield */
	    fdpad(npml,pad1,pad2, srce[is],Bx[its],Bz[its]);	    

	    (void) umfpack_zl_solve (UMFPACK_A, 
				     NULL, NULL, NULL, NULL, 
				     Xx[its], Xz[its], Bx[its], Bz[its], 
				     Numeric[its], Control, NULL);
	    
	    fdcut(npml,pad1,pad2, srce[is],Xx[its],Xz[its]);

	    /* receiver wavefield */
	    fdpad(npml,pad1,pad2, recv[is],Bx[its],Bz[its]);	    
	    
	    (void) umfpack_zl_solve (UMFPACK_At, 
				     NULL, NULL, NULL, NULL, 
				     Xx[its], Xz[its], Bx[its], Bz[its], 
				     Numeric[its], Control, NULL);
	    	    
	    fdcut(npml,pad1,pad2, recv[is],Xx[its],Xz[its]);

	    /* imaging condition */
	    for (ih=-nh; ih < nh+1; ih++) {
		for (j=0; j < n2; j++) {
		    for (i=0; i < n1; i++) {
			if (j-abs(ih) >= 0 && j+abs(ih) < n2) {
			    image[its][ih+nh][j][i] += crealf(conjf(srce[is][j-ih][i])*recv[is][j+ih][i]);
			}
		    }
		}
	    }
	}

	(void) umfpack_zl_free_symbolic (&Symbolic);
	for (its=0; its < uts; its++) {
	    (void) umfpack_zl_free_numeric (&Numeric[its]);
	}

	free(Ti); free(Tj); free(Tx); free(Tz);
	free(Ap); free(Ai); free(Map);
	free(Ax); free(Az);
	for (its=0; its < uts; its++) {
	    free(Bx[its]); free(Bz[its]); free(Xx[its]); free(Xz[its]);
	}

	/* write wavefields to temporary file */
	sf_fslice_put(sfile,iw,srce[0][0]);
	sf_fslice_put(rfile,iw,recv[0][0]);
    }

    /* rewind file stream */
    rewind(sf_filestream(source));
    rewind(sf_filestream(data));

    /* output extended image */
    for (ih=-nh; ih < nh+1; ih++) {
	for (j=0; j < n2; j++) {
	    for (i=0; i < n1; i++) {
		image0[(ih+nh)*n1*n2+j*n1+i] = image[0][ih+nh][j][i];
	    }
	}
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(uts) private(j,i,its)
    for (ih=-nh; ih < nh+1; ih++) {
	for (j=0; j < n2; j++) {
	    for (i=0; i < n1; i++) {
		for (its=1; its < uts; its++) {
		    image0[(ih+nh)*n1*n2+j*n1+i] += image[its][ih+nh][j][i];
		}
	    }
	}
    }
#endif
}
