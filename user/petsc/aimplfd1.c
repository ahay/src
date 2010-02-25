/* Implicit solver for 1-D acoustic wave equation with PETSc. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include <petscksp.h>
/*^*/
#include <rsf.h>
/*^*/

#include "aimplfd1.h"


#ifndef _aimplfd1_h

typedef struct PETScAimplFD1 *sf_petsc_aimplfd1;
/* abstract data type */
/*^*/
#endif

struct PETScAimplFD1 {
    PetscInt N;
    MPI_Comm comm;
    KSP Solver;
    Mat A;
    Vec Rhs, Ut, Ut1, Ut2;
    float dx, dt;
};
/* concrete data type */

#define CHKERR CHKERRABORT(aimplfd->comm,ierr)

sf_petsc_aimplfd1 sf_petsc_aimplfd1_init (int nx, float dx, float dt,
                                          float *v, int niter)
/*< Initialize implicit F-D object. >*/
{
    PetscInt J, K;
    int i;
    float ro, rtx2;
    PetscErrorCode ierr;
    PetscScalar Ro12, Rom, Val;

    sf_petsc_aimplfd1 aimplfd;

    aimplfd = (sf_petsc_aimplfd1)sf_alloc (1, sizeof (struct PETScAimplFD1));

    aimplfd->N = nx;
    aimplfd->dx = dx;
    aimplfd->dt = dt;
    aimplfd->comm = MPI_COMM_WORLD;

    /* Matrix */
    ierr = MatCreate (aimplfd->comm, &aimplfd->A); CHKERR;
    /* Matrix width/heigh */
    ierr = MatSetSizes (aimplfd->A, PETSC_DECIDE, PETSC_DECIDE,
                        aimplfd->N, aimplfd->N); CHKERR;
    ierr = MatSetType (aimplfd->A, MATMPIAIJ); CHKERR;
    /* Solver type */
    ierr = KSPCreate (aimplfd->comm, &aimplfd->Solver); CHKERR;
    ierr = KSPSetOperators (aimplfd->Solver, aimplfd->A, aimplfd->A, SAME_NONZERO_PATTERN); CHKERR;
    ierr = KSPSetType (aimplfd->Solver, KSPGMRES); CHKERR;
    ierr = KSPSetTolerances (aimplfd->Solver, PETSC_DEFAULT, PETSC_DEFAULT,
                             PETSC_DEFAULT, niter); CHKERR;
    /* Solution and r.h.s. vectors */
    ierr = VecCreateMPI (aimplfd->comm, PETSC_DECIDE, aimplfd->N, &aimplfd->Rhs); CHKERR;
    Val = 0.0;
    ierr = VecSet (aimplfd->Rhs, Val); CHKERR;
    ierr = VecDuplicate (aimplfd->Rhs, &aimplfd->Ut); CHKERR;
    ierr = VecDuplicate (aimplfd->Rhs, &aimplfd->Ut1); CHKERR;
    ierr = VecDuplicate (aimplfd->Rhs, &aimplfd->Ut2); CHKERR;

    /* Build matrix */
    ierr = MatAssemblyBegin (aimplfd->A, MAT_FINAL_ASSEMBLY); CHKERR;
    rtx2 = (dt*dt)/(dx*dx);
    /* Left B.C. */
    J = 0;
    K = 0;
    ro = v[0]*rtx2;
    Ro12 = 1 + 2*ro;
    Rom = -ro;
    ierr = MatSetValues (aimplfd->A, 1, &J, 1, &K, &Ro12, INSERT_VALUES); CHKERR;
    K = K + 1;
    ierr = MatSetValues (aimplfd->A, 1, &J, 1, &K, &Rom, INSERT_VALUES); CHKERR;
    for (i = 1; i < (nx - 1); i++) {
        J = i;
        K = i - 1;
        ro = v[i]*rtx2;
        /* Coefficients */
        Ro12 = 1 + 2*ro;
        Rom = -ro;
        ierr = MatSetValues (aimplfd->A, 1, &J, 1, &K, &Rom, INSERT_VALUES); CHKERR;
        K = K + 1;
        ierr = MatSetValues (aimplfd->A, 1, &J, 1, &K, &Ro12, INSERT_VALUES); CHKERR;
        K = K + 1;
        ierr = MatSetValues (aimplfd->A, 1, &J, 1, &K, &Rom, INSERT_VALUES); CHKERR;
    }
    /* Right B.C. */
    J = nx - 1;
    K = nx - 2;
    ro = v[nx - 1]*rtx2;
    Ro12 = 1 + 2*ro;
    Rom = -ro;
    ierr = MatSetValues (aimplfd->A, 1, &J, 1, &K, &Rom, INSERT_VALUES); CHKERR;
    K = K + 1;
    ierr = MatSetValues (aimplfd->A, 1, &J, 1, &K, &Ro12, INSERT_VALUES); CHKERR;
    /* Finish matrix assembly */
    ierr = MatAssemblyEnd (aimplfd->A, MAT_FINAL_ASSEMBLY); CHKERR;
    ierr = MatSetOption (aimplfd->A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE); CHKERR;

    return aimplfd;
}

void sf_petsc_aimplfd1_next_step (sf_petsc_aimplfd1 aimplfd)
/*< Calculate next time step. >*/
{
    PetscErrorCode ierr;
    PetscScalar Val;

    /* R.H.S. */
    Val = -1.0;
    VecScale (aimplfd->Ut2, Val);
    Val = 2.0;
    VecWAXPY (aimplfd->Rhs, Val, aimplfd->Ut1, aimplfd->Ut2);

    /* Solve the system */
    ierr = KSPSolve (aimplfd->Solver, aimplfd->Rhs, aimplfd->Ut); CHKERR;

    /* Move time steps */
    ierr = VecCopy (aimplfd->Ut1, aimplfd->Ut2); CHKERR;
    ierr = VecCopy (aimplfd->Ut, aimplfd->Ut1); CHKERR;
}

void sf_petsc_aimplfd1_add_source (sf_petsc_aimplfd1 aimplfd, float f, int ix)
/*< Inject source at a specific location. >*/
{
    PetscInt J;
    PetscErrorCode ierr;
    PetscScalar Val;

    ierr = VecAssemblyBegin (aimplfd->Ut1); CHKERR;
    Val = f*0.5*aimplfd->dt*aimplfd->dt;
    J = ix;
    ierr = VecSetValue (aimplfd->Ut1, J, Val, ADD_VALUES); CHKERR;
    ierr = VecAssemblyEnd (aimplfd->Ut1); CHKERR;
}

void sf_petsc_aimplfd1_get_wavefield (sf_petsc_aimplfd1 aimplfd, float *u)
/*< Get wavefield values at the current time step. >*/
{
    int i, nx;
    PetscInt J;
    PetscErrorCode ierr;
    PetscScalar Val;
    Vec Uout;
    VecScatter Uctx;

    ierr = VecScatterCreateToZero (aimplfd->Ut, &Uctx, &Uout); CHKERR;
    ierr = VecScatterBegin (Uctx, aimplfd->Ut, Uout, INSERT_VALUES, SCATTER_FORWARD); CHKERR;
    nx = aimplfd->N;
    for (i = 0; i < nx; i++) {
        J = i;
        ierr = VecGetValues (Uout, 1, &J, &Val); CHKERR;
        u[i] = Val;
    }
    ierr = VecScatterEnd (Uctx, aimplfd->Ut, Uout, INSERT_VALUES, SCATTER_FORWARD); CHKERR;
    ierr = VecScatterDestroy (Uctx); CHKERR;
    ierr = VecDestroy (Uout); CHKERR;
}

void sf_petsc_aimplfd1_destroy (sf_petsc_aimplfd1 aimplfd)
/*< Destroy the object. >*/
{
    PetscErrorCode ierr;

    ierr = KSPDestroy (aimplfd->Solver); CHKERR;
    ierr = VecDestroy (aimplfd->Rhs); CHKERR;
    ierr = VecDestroy (aimplfd->Ut); CHKERR;
    ierr = VecDestroy (aimplfd->Ut1); CHKERR;
    ierr = VecDestroy (aimplfd->Ut2); CHKERR;

    free (aimplfd);
}

