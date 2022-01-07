/* Interface for LU factorization */
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


/*^*/
#ifndef SuiteSparse_long
#ifdef UF_long
#define SuiteSparse_long UF_long
#endif
#endif
/*^*/


#include "fdprep.h"
#include "sparsesolver.h"


static double **Bx, **Bz, **Xx, **Xz;
static void **Numeric, *Symbolic;
static double Control[UMFPACK_CONTROL];
static double *Tx, *Tz, *Ax, *Az;
static SuiteSparse_long *Ti, *Tj, *Ap, *Ai, *Map;
static SuiteSparse_long n, nz;

void sparse_init(int uts, int pad1, int pad2)
/*< sparse init >*/
{
    int its; 

    n  = fdprep_n (pad1,pad2);
	nz = fdprep_nz(pad1,pad2);

    Numeric = (void**) sf_alloc(uts,sizeof(void*)); 

    umfpack_zl_defaults(Control);
    Control [UMFPACK_IRSTEP]=0; 

    Ti = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long));
    Tj = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long)); 
    Ap = (SuiteSparse_long*) sf_alloc(n+1,sizeof(SuiteSparse_long));
    Ai = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long));
    Map = (SuiteSparse_long*) sf_alloc(nz,sizeof(SuiteSparse_long)); 

    Ax = (double*) sf_alloc(nz,sizeof(double)); 
    Az = (double*) sf_alloc(nz,sizeof(double));
    Tx = (double*) sf_alloc(nz,sizeof(double));
    Tz = (double*) sf_alloc(nz,sizeof(double)); 

    Bx = (double**) sf_alloc(uts,sizeof(double*));
    Bz = (double**) sf_alloc(uts,sizeof(double*));
    Xx = (double**) sf_alloc(uts,sizeof(double*));
    Xz = (double**) sf_alloc(uts,sizeof(double*));

    for (its=0; its<uts; its++ ) { 
        Bx[its] = (double*) sf_alloc(n,sizeof(double));
        Bz[its] = (double*) sf_alloc(n,sizeof(double));
        Xx[its] = (double*) sf_alloc(n,sizeof(double));
        Xz[its] = (double*) sf_alloc(n,sizeof(double));
    }
}

void sparse_factor(double omega, int n1, int n2, float d1, float d2,
                float **v, int npml, int pad1, int pad2)
/*< sparse factor>*/
{
    
    fdprep(omega, n1, n2, d1, d2, v,
            npml, pad1, pad2, Ti, Tj, Tx, Tz);


    umfpack_zl_triplet_to_col(n, n, nz, 
        Ti, Tj, Tx, Tz, Ap, Ai, Ax, Az, Map);


    umfpack_zl_symbolic(n, n, Ap, Ai, Ax, Az, 
        &Symbolic, Control, NULL);

    
    umfpack_zl_numeric(Ap, Ai, Ax, Az, Symbolic, &Numeric[0],
        Control, NULL);

}

void sparse_solve(int npml, int pad1, int pad2, sf_complex ***f, bool hermite, int ns, int uts)
/*< sparse solver>*/
{
    int is, its;

    /* loop over shots */
//#ifdef _OPENMP
//#pragma omp parallel for num_threads(uts) private(its)
//#endif
    for ( is=0; is < ns; is++ ) {
//#ifdef _OPENMP
//	    its = omp_get_thread_num();
//#else
	    its = 0;
//#endif

        fdpad(npml, pad1, pad2, f[is], Bx[its], Bz[its]);

        umfpack_zl_solve(hermite?UMFPACK_At: UMFPACK_A, 
            NULL, NULL, NULL, NULL,
            Xx[its], Xz[its], Bx[its], Bz[its],
            Numeric[its], Control, NULL);
        fdcut(npml, pad1, pad2, f[is], Xx[its], Xz[its]);
    }
}


void sparse_free(int uts)
/*< sparse free>*/
{
    int its;

    (void) umfpack_zl_free_symbolic(&Symbolic);
    for (its=0; its<uts; its++) {
        (void) umfpack_zl_free_numeric(&Numeric[its]);
    }

    free(Ti); free(Tj); free(Tx); free(Tz);
    free(Ap); free(Ai); free(Map);
    free(Ax); free(Az);

    for (its=0; its< uts; its++ ) {
        free(Bx[its]); free(Bz[its]); free(Xx[its]); free(Xz[its]);
    }
}

