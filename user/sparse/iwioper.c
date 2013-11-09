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

static float **vel, d1, d2, ow, dw;
static float ***wght, **prec;
static float **tempx, **tempr;
static int n1, n2, nh, ns, nw, npml;
static sf_file sfile, rfile;
static SuiteSparse_long *Ti, *Tj;
static SuiteSparse_long *Ap, *Ai, *Map;
static void *Symbolic, **Numeric;
static double Control[UMFPACK_CONTROL];
static double *Tx, *Tz;
static double *Ax, *Az, **Xx, **Xz, **Bx, **Bz;
static sf_complex ***us, ***ur, ***as, ***ar;
static bool load, mass;
static int uts, ss[3];
static char *datapath, *insert, *append;
static size_t srclen, inslen;
static char *order;

void adjsrce(sf_complex **recv /* receiver wavefield */,
	     sf_complex **adjs /* adjoint-source */,
	     float *dm, float *di, bool adj)
/* assemble ajoint-source */
{
    int i, j, ih;

    if (adj) {
	for (ih=-nh; ih < nh+1; ih++) {
	    for (j=0; j < n2; j++) {
		for (i=0; i < n1; i++) {
		    if (j+2*ih >= 0 && j+2*ih < n2) {
			adjs[j][i] += recv[j+2*ih][i]
			    *(wght==NULL? 1.: wght[ih+nh][j+ih][i])
			    *di[(ih+nh)*ss[2]+(j+ih)*ss[1]+i];
		    }
		}
	    }
	}
    } else {
	for (j=0; j < n2; j++) {
	    for (i=0; i < n1; i++) {
		adjs[j][i] = recv[j][i]
		    *(prec==NULL? 1.: prec[j][i])
		    *dm[j*ss[1]+i];
	    }
	}
    }
}

void adjrecv(sf_complex **srce /* source wavefield */,
	     sf_complex **adjr /* adjoint-receiver */,
	     float *dm, float *di, bool adj)
/* assemble ajoint-receiver */
{
    int i, j, ih;

    if (adj) {
	for (ih=-nh; ih < nh+1; ih++) {
	    for (j=0; j < n2; j++) {
		for (i=0; i < n1; i++) {
		    if (j-2*ih >= 0 && j-2*ih < n2) {
			adjr[j][i] += srce[j-2*ih][i]
			    *(wght==NULL? 1.: wght[ih+nh][j-ih][i])
			    *di[(ih+nh)*ss[2]+(j-ih)*ss[1]+i];
		    }
		}
	    }
	}
    } else {
	for (j=0; j < n2; j++) {
	    for (i=0; i < n1; i++) {
		adjr[j][i] = srce[j][i]
		    *(prec==NULL? 1.: prec[j][i])
		    *dm[j*ss[1]+i];
	    }
	}
    }
}

void adjclean(sf_complex **adjs,
	      sf_complex **adjr)
/* clean-up */
{
    int i, j;

    for (j=0; j < n2; j++) {
	for (i=0; i < n1; i++) {
	    adjs[j][i] = sf_cmplx(0.,0.);
	    adjr[j][i] = sf_cmplx(0.,0.);
	}
    }
}

sf_complex iwiavrg(sf_complex **wave,
		   int i, int j)
/* mass term */
{
    sf_complex temp;

    switch (order[0]) {
	case '5':
	    return wave[j][i];
	    break;

	case '9':
	    return wave[j][i];
	    break;

	case 'j':
	    temp = 0.6296*wave[j][i];
	    if (i > 0               ) temp += 0.0942*wave[j][i-1];
	    if (i < n1-1            ) temp += 0.0942*wave[j][i+1];
	    if (            j > 0   ) temp += 0.0942*wave[j-1][i];
	    if (            j < n2-1) temp += 0.0942*wave[j+1][i];
	    if (i > 0    && j > 0   ) temp -= 0.0016*wave[j-1][i-1];
	    if (i > 0    && j < n2-1) temp -= 0.0016*wave[j+1][i-1];
	    if (i < n1-1 && j < n2-1) temp -= 0.0016*wave[j+1][i+1];
	    if (i < n1-1 && j > 0   ) temp -= 0.0016*wave[j-1][i+1];
	    return temp;
	    break;

	case 'c':
	    temp = 0.363276*wave[j][i];
	    if (i > 0               ) temp += 0.108598*wave[j][i-1];
	    if (i < n1-1            ) temp += 0.108598*wave[j][i+1];
	    if (            j > 0   ) temp += 0.108598*wave[j-1][i];
	    if (            j < n2-1) temp += 0.108598*wave[j+1][i];
	    if (i > 0    && j > 0   ) temp += 0.0424801*wave[j-1][i-1];
	    if (i > 0    && j < n2-1) temp += 0.0424801*wave[j+1][i-1];
	    if (i < n1-1 && j < n2-1) temp += 0.0424801*wave[j+1][i+1];
	    if (i < n1-1 && j > 0   ) temp += 0.0424801*wave[j-1][i+1];
	    if (i > 1               ) temp += 0.0041487*wave[j][i-2];
	    if (i < n1-2            ) temp += 0.0041487*wave[j][i+2];
	    if (            j > 1   ) temp += 0.0041487*wave[j-2][i];
	    if (            j < n2-2) temp += 0.0041487*wave[j+2][i];
	    if (i > 1    && j > 1   ) temp += 0.000206312*wave[j-2][i-2];
	    if (i > 1    && j < n2-2) temp += 0.000206312*wave[j+2][i-2];
	    if (i < n1-2 && j < n2-2) temp += 0.000206312*wave[j+2][i+2];
	    if (i < n1-2 && j > 1   ) temp += 0.000206312*wave[j-2][i+2];
	    if (i > 1    && j > 0   ) temp += 0.00187765*wave[j-1][i-2];
	    if (i > 0    && j < n2-2) temp += 0.00187765*wave[j+2][i-1];
	    if (i < n1-2 && j < n2-1) temp += 0.00187765*wave[j+1][i+2];
	    if (i < n1-1 && j > 1   ) temp += 0.00187765*wave[j-2][i+1];
	    if (i > 1    && j < n2-1) temp += 0.00188342*wave[j+1][i-2];
	    if (i < n1-1 && j < n2-2) temp += 0.00188342*wave[j+2][i+1];
	    if (i < n1-2 && j > 0   ) temp += 0.00188342*wave[j-1][i+2];
	    if (i > 0    && j > 1   ) temp += 0.00188342*wave[j-2][i-1];
	    return temp;
	    break;

	default:
	    sf_error("Fail to load discretization scheme.");
    }

    return sf_cmplx(0.,0.);
}

void iwiadd(double omega,	     
	    sf_complex **srce /* source */,
	    sf_complex **recv /* receiver */,
	    sf_complex **adjs /* adjoint-source */,
	    sf_complex **adjr /* adjoint-receiver */,
	    float *dm, float *di, bool adj)
/* assemble */
{
    int i, j, ih;

    if (adj) {
	for (j=0; j < n2; j++) {
	    for (i=0; i < n1; i++) {    
		if (mass) {
		    dm[j*ss[1]+i] -= omega*omega
			*(prec==NULL? 1.: prec[j][i])*crealf(
			    conjf(iwiavrg(srce,i,j))*adjs[j][i]+
			    iwiavrg(recv,i,j)*conjf(adjr[j][i]));
		} else {
		   dm[j*ss[1]+i] -= omega*omega
		       *(prec==NULL? 1.: prec[j][i])*crealf(
			   conjf(srce[j][i])*adjs[j][i]+
			   recv[j][i]*conjf(adjr[j][i]));
		} 
	    }
	}
    } else {
	for (ih=-nh; ih < nh+1; ih++) {
	    for (j=0; j < n2; j++) {
		for (i=0; i < n1; i++) {
		    if (j-abs(ih) >= 0 && j+abs(ih) < n2) {
			if (mass) {
			    di[(ih+nh)*ss[2]+j*ss[1]+i] -= omega*omega
				*(wght==NULL? 1.: wght[ih+nh][j][i])*crealf(
				    recv[j+ih][i]*conj(iwiavrg(adjr,i,j-ih))+
				    conjf(srce[j-ih][i])*iwiavrg(adjs,i,j+ih));
			} else {
			    di[(ih+nh)*ss[2]+j*ss[1]+i] -= omega*omega
				*(wght==NULL? 1.: wght[ih+nh][j][i])*crealf(
				    recv[j+ih][i]*conj(adjr[j-ih][i])+
				    conjf(srce[j-ih][i])*adjs[j+ih][i]);
			}
		    }
		}
	    }
	}
    }
}

void iwi_init(int npml0, 
	      int nn1, int nn2, 
	      float dd1, float dd2,
	      int nh0, int ns0, 
	      float ow0, float dw0, int nw0,
	      sf_file us0, sf_file ur0,
	      bool load0, char *datapath0,
	      int uts0,
	      char *order0, bool mass0)
/*< initialize >*/
{
    npml = npml0;
    
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
    
    uts = uts0;

    order = order0;
    fdprep_order(order);

    mass = mass0;

    ss[0] = 1; ss[1] = n1; ss[2] = n1*n2;

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
    us = sf_complexalloc3(n1,n2,ns);
    ur = sf_complexalloc3(n1,n2,ns);
    as = sf_complexalloc3(n1,n2,ns);
    ar = sf_complexalloc3(n1,n2,ns);

    if (load) {
	Ti = NULL; Tj = NULL; Tx = NULL; Tz = NULL; 
	Ap = NULL; Ai = NULL; Map = NULL; Ax = NULL; Az = NULL;
    }

    Bx = (double**) sf_alloc(uts,sizeof(double*));
    Bz = (double**) sf_alloc(uts,sizeof(double*));
    Xx = (double**) sf_alloc(uts,sizeof(double*));
    Xz = (double**) sf_alloc(uts,sizeof(double*));

    Numeric = (void**) sf_alloc(uts,sizeof(void*));

    tempx = sf_floatalloc2(n1*n2,uts);
    tempr = sf_floatalloc2(n1*n2*(2*nh+1),uts);

    /* LU control */
    umfpack_zl_defaults (Control);
    Control [UMFPACK_IRSTEP] = 0;
}

void iwi_free()
/*< free memories >*/
{
    int its;

    for (its=0; its < uts; its++) {
	free(tempx[its]); free(tempr[its]);
    }
    free(tempx); free(tempr);
}

void iwi_set(float **vel0,
	     float ***wght0,
	     float **prec0)
/*< set velocity, weight and preconditioner >*/
{
    vel = vel0;
    wght = wght0;
    prec = prec0;
}

void iwi_oper(bool adj, bool add, int nx, int nr, float *x, float *r)
/*< linear operator >*/
{
    int iw, is, its, i;
    int pad1, pad2;
    SuiteSparse_long n, nz;
    double omega;

    sf_adjnull(adj,add,nx,nr,x,r);

    /* PML padding */
    pad1 = n1+2*npml;
    pad2 = n2+2*npml;

    n  = fdprep_n (pad1,pad2);
    nz = fdprep_nz(pad1,pad2);

    /* loop over frequency */
    for (iw=0; iw < nw; iw++) {
	omega = (double) 2.*SF_PI*(ow+iw*dw);	

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

	/* LU file (append _lu* after velocity file) */
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

#ifdef _OPENMP
	    (void) umfpack_zl_save_numeric (Numeric[0], append);
	    
	    for (its=1; its < uts; its++) {
		(void) umfpack_zl_load_numeric (&Numeric[its], append);
	    }
	    
	    (void) remove (append);
	    (void) remove ("numeric.umf");
#endif
	} else {
	    /* load Numeric */
	    for (its=0; its < uts; its++) {
		(void) umfpack_zl_load_numeric (&Numeric[its], append);
	    }
	}

	if (load) free(append);

	/* background wavefields */
	sf_complexread(us[0][0],n1*n2*ns,sfile);	    
	sf_complexread(ur[0][0],n1*n2*ns,rfile);

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

	    /* adjoint source */
	    adjsrce(ur[is],as[is], x,r,adj);

	    fdpad(npml,pad1,pad2, as[is],Bx[its],Bz[its]);
	    
	    (void) umfpack_zl_solve (UMFPACK_At, 
				     NULL, NULL, NULL, NULL, 
				     Xx[its], Xz[its], Bx[its], Bz[its], 
				     Numeric[its], Control, NULL);
	    	    
	    fdcut(npml,pad1,pad2, as[is],Xx[its],Xz[its]);

	    /* adjoint receiver */
	    adjrecv(us[is],ar[is], x,r,adj);

	    fdpad(npml,pad1,pad2, ar[is],Bx[its],Bz[its]);	    

	    (void) umfpack_zl_solve (UMFPACK_A, 
				     NULL, NULL, NULL, NULL, 
				     Xx[its], Xz[its], Bx[its], Bz[its], 
				     Numeric[its], Control, NULL);
	    
	    fdcut(npml,pad1,pad2, ar[is],Xx[its],Xz[its]);

	    /* assemble */
	    iwiadd(omega, us[is],ur[is],as[is],ar[is], tempx[its],tempr[its],adj);

	    /* clean up */
	    if (adj) adjclean(as[is],ar[is]);
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

    /* rewind file stream */
    rewind(sf_filestream(sfile));
    rewind(sf_filestream(rfile));

#ifdef _OPENMP
    if (adj) {
#pragma omp parallel for num_threads(uts) private(its)
	for (i=0; i < n1*n2; i++) {
	    for (its=0; its < uts; its++) {
		x[i] += tempx[its][i];
		tempx[its][i] = 0.;
	    }
	}
    } else {
#pragma omp parallel for num_threads(uts) private(its)
	for (i=0; i < n1*n2*(2*nh+1); i++) {
	    for (its=0; its < uts; its++) {
		r[i] += tempr[its][i];
		tempr[its][i] = 0.;
	    }
	}
    }
#else
    if (adj) {
	for (i=0; i < n1*n2; i++) {
	    x[i] = tempx[0][i];
	    tempx[0][i] = 0.;
	}
    } else {
	for (i=0; i < n1*n2*(2*nh+1); i++) {
	    r[i] = tempr[0][i];
	    tempr[0][i] = 0.;
	}
    }
#endif
}
