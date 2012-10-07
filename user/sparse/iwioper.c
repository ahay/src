/* Image-domain waveform tomography (linear operator). */
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
#include "iwioper.h"

static float eps, **vel, d1, d2, ow, dw, ***wght;
static int n1, n2, npml, pad1, pad2, nh, ns, nw, ss[3];
static sf_file sfile, rfile;
static SuiteSparse_long n, nz, *Ti, *Tj;
static SuiteSparse_long status, *Ap, *Ai, *Map;
static void *Symbolic, *Numeric;
static double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
static double *Tx, *Tz;
static double *Ax, *Az, *Xx, *Xz, *Bx, *Bz;
static sf_complex ***us, ***ur, ***as, ***ar;
static bool load;
static char *datapath, *insert, *append;
static size_t srclen, inslen;

void adjsrce(sf_complex ***recv /* receiver wavefield */,
	     float *rhs /* right-hand side */,
	     sf_complex ***adjs /* adjoint-source */,
	     bool adj)
/* assemble ajoint-source */
{
    int is, i, j, ih;

    if (adj) {
#ifdef _OPENMP    
#pragma omp parallel for private(ih,j,i)
#endif
	for (is=0; is < ns; is++) {
	    for (ih=-nh; ih < nh+1; ih++) {
		for (j=0; j < n2; j++) {
		    for (i=0; i < n1; i++) {
			if (j+2*ih >= 0 && j+2*ih < n2) {
			    adjs[is][j][i] += recv[is][j+2*ih][i]
				*(wght==NULL? 1.: wght[ih+nh][j+ih][i])
				*rhs[(ih+nh)*ss[2]+(j+ih)*ss[1]+i];
			}
		    }
		}
	    }
	}
    } else {
#ifdef _OPENMP    
#pragma omp parallel for private(j,i)
#endif	
	for (is=0; is < ns; is++) {
	    for (j=0; j < n2; j++) {
		for (i=0; i < n1; i++) {
		    adjs[is][j][i] = recv[is][j][i]*rhs[j*ss[1]+i];
		}
	    }
	}
    }
}

void adjrecv(sf_complex ***srce /* source wavefield */,
	     float *rhs /* right-hand side */,
	     sf_complex ***adjr /* adjoint-receiver */,
	     bool adj)
/* assemble ajoint-receiver */
{
    int is, i, j, ih;

    if (adj) {
#ifdef _OPENMP    
#pragma omp parallel for private(ih,j,i)
#endif
	for (is=0; is < ns; is++) {
	    for (ih=-nh; ih < nh+1; ih++) {
		for (j=0; j < n2; j++) {
		    for (i=0; i < n1; i++) {
			if (j-2*ih >= 0 && j-2*ih < n2) {
			    adjr[is][j][i] += srce[is][j-2*ih][i]
				*(wght==NULL? 1.: wght[ih+nh][j-ih][i])
				*rhs[(ih+nh)*ss[2]+(j-ih)*ss[1]+i];
			}
		    }
		}
	    }
	}
    } else {
#ifdef _OPENMP    
#pragma omp parallel for private(j,i)
#endif
	for (is=0; is < ns; is++) {
	    for (j=0; j < n2; j++) {
		for (i=0; i < n1; i++) {
		    adjr[is][j][i] = srce[is][j][i]*rhs[j*ss[1]+i];
		}
	    }
	}
    }
}

void adjclean(sf_complex ***adjs,
	      sf_complex ***adjr)
/* clean-up */
{
    int is, i, j;

#ifdef _OPENMP    
#pragma omp parallel for private(j,i)
#endif
    for (is=0; is < ns; is++) {
	for (j=0; j < n2; j++) {
	    for (i=0; i < n1; i++) {
		adjs[is][j][i] = 0.;
		adjr[is][j][i] = 0.;
	    }
	}
    }
}

void iwiadd(double omega,	     
	    sf_complex ***srce /* source */,
	    sf_complex ***recv /* receiver */,
	    sf_complex ***adjs /* adjoint-source */,
	    sf_complex ***adjr /* adjoint-receiver */,
	    float *dm, float *di, bool adj)
/* assemble */
{
    int is, i, j, ih;

    if (adj) {
	for (is=0; is < ns; is++) {
	    for (j=0; j < n2; j++) {
		for (i=0; i < n1; i++) {    
		    dm[j*ss[1]+i] -= omega*omega*creal(
			conjf(srce[is][j][i])*adjs[is][j][i]+
			recv[is][j][i]*conjf(adjr[is][j][i]));
		}
	    }
	}
    } else {
#ifdef _OPENMP    
#pragma omp parallel for private(is,j,i)
#endif
	for (ih=-nh; ih < nh+1; ih++) {
	    for (is=0; is < ns; is++) {
		for (j=0; j < n2; j++) {
		    for (i=0; i < n1; i++) {
			if (j-abs(ih) >= 0 && j+abs(ih) < n2) {
			    di[(ih+nh)*ss[2]+j*ss[1]+i] -= omega*omega
				*(wght==NULL? 1.: wght[ih+nh][j][i])*creal(
				    recv[is][j+ih][i]*conj(adjr[is][j-ih][i])+
				    conjf(srce[is][j-ih][i])*adjs[is][j+ih][i]);
			}
		    }
		}
	    }
	}
    }
}

void iwi_init(int npw, float eps0,
	      int nn1, int nn2, 
	      float dd1, float dd2,
	      int nh0, int ns0, 
	      float ow0, float dw0, int nw0,
	      sf_file us0, sf_file ur0,
	      bool load0, char *datapath0)
/*< initialize >*/
{
    eps = eps0;

    n1 = nn1;
    n2 = nn2;
    d1 = dd1;
    d2 = dd2;

    nh = nh0; ns = ns0;
    ow = ow0; dw = dw0; nw = nw0;

    sfile = us0;
    rfile = ur0;

    load = load0;
    datapath = datapath0;

    ss[0] = 1; ss[1] = n1; ss[2] = n1*n2;

    /* LU file */
    if (load) {
	srclen = strlen(datapath);
	insert = sf_charalloc(6);
    } else {
	insert = NULL;
	append = NULL;
    }

    /* prepare PML and LU */
    npml = npw*2;
    pad1 = n1+2*npml;
    pad2 = n2+2*npml;

    n = (pad1-2)*(pad2-2);
    nz = 5*(pad1-2)*(pad2-2)-2*(pad1-4)-2*(pad2-4)-8;

    /* allocate temporary space */
    us = sf_complexalloc3(n1,n2,ns);
    ur = sf_complexalloc3(n1,n2,ns);
    as = sf_complexalloc3(n1,n2,ns);
    ar = sf_complexalloc3(n1,n2,ns);

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
    } else {
	Ti = NULL; Tj = NULL; Tx = NULL; Tz = NULL; 
	Ap = NULL; Ai = NULL; Map = NULL; Ax = NULL; Az = NULL;
    }

    Bx = (double*) sf_alloc(n,sizeof(double));
    Bz = (double*) sf_alloc(n,sizeof(double));
    Xx = (double*) sf_alloc(n,sizeof(double));
    Xz = (double*) sf_alloc(n,sizeof(double));

    /* turn off iterative refinement */
    umfpack_zl_defaults (Control);
    Control [UMFPACK_IRSTEP] = 0;
}

void iwi_set(float **vel0,
	     float ***wght0)
/*< set velocity and weight >*/
{
    vel = vel0;
    wght = wght0;
}

void iwi_oper(bool adj, bool add, int nx, int nr, float *x, float *r)
/*< linear operator >*/
{
    int iw, is;
    double omega;

    sf_adjnull(adj,add,nx,nr,x,r);

    /* loop over frequency */
    for (iw=0; iw < nw; iw++) {
	omega = (double) 2.*SF_PI*(ow+iw*dw);

	/* load LU */
	if (load) {
	    sprintf(insert,"_lu%d",iw);
	    inslen = strlen(insert);
	    
	    append = malloc(srclen+inslen+1);
	    
	    memcpy(append,datapath,srclen-5);
	    memcpy(append+srclen-5,insert,inslen);
	    memcpy(append+srclen-5+inslen,datapath+srclen-5,5+1);
	}

	if (!load) {
	    /* assemble matrix */
	    fdprep(omega, eps, 
		   n1, n2, d1, d2, vel,
		   npml, pad1, pad2, n, nz, 
		   Ti, Tj, Tx, Tz);
	    
	    status = umfpack_zl_triplet_to_col (n, n, nz, 
						Ti, Tj, Tx, Tz, 
						Ap, Ai, Ax, Az, Map);
	    
	    /* LU */
	    status = umfpack_zl_symbolic (n, n, 
					  Ap, Ai, Ax, Az, 
					  &Symbolic, Control, Info);
	    
	    status = umfpack_zl_numeric (Ap, Ai, Ax, Az, 
					 Symbolic, &Numeric, 
					 Control, Info);
	} else {
	    /* load Numeric */
	    status = umfpack_zl_load_numeric (&Numeric, append);
	}

	if (load) free(append);

	/* background wavefields */
	sf_complexread(us[0][0],n1*n2*ns,sfile);	    
	sf_complexread(ur[0][0],n1*n2*ns,rfile);

	/* adjoint wavefields */
	adjsrce(ur, r,as, adj);
	adjrecv(us, r,ar, adj);

	for (is=0; is < ns; is++) {
	    fdpad(npml,pad1,pad2, as[is],Bx,Bz);

	    status = umfpack_zl_solve (UMFPACK_At, 
				       NULL, NULL, NULL, NULL, 
				       Xx, Xz, Bx, Bz, 
				       Numeric, Control, Info);

	    fdcut(npml,pad1,pad2, as[is],Xx,Xz);

	    fdpad(npml,pad1,pad2, ar[is],Bx,Bz);

	    status = umfpack_zl_solve (UMFPACK_A, 
				       NULL, NULL, NULL, NULL, 
				       Xx, Xz, Bx, Bz, 
				       Numeric, Control, Info);

	    fdcut(npml,pad1,pad2, ar[is],Xx,Xz);
	}

	/* assemble */
	iwiadd(omega, us,ur,as,ar, x,r,adj);
	
	/* clean up */
	if (adj) adjclean(as,ar);

	if (!load) (void) umfpack_zl_free_symbolic (&Symbolic);
	(void) umfpack_zl_free_numeric (&Numeric);
    }
}
