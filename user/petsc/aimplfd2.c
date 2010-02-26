/* Implicit solver for 2-D acoustic wave equation with PETSc. */
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

#include "aimplfd2.h"


#ifndef _aimplfd2_h

typedef struct PETScAimplFD2 *sf_petsc_aimplfd2;
/* abstract data type */
/*^*/
#endif

struct PETScAimplFD2 {
    PetscInt Nx, Nz;
    MPI_Comm comm;
    KSP Solver;
    Mat A;
    Vec Rhs, Ut, Ut1, Ut2;
    float dz, dx, dt;
};
/* concrete data type */

#define CHKERR CHKERRABORT(aimplfd->comm,ierr)

static void sf_petsc_aimplfd2_set_mat (sf_petsc_aimplfd2 aimplfd,
                                       PetscInt J, PetscInt K, PetscInt N, PetscScalar C)
{
    PetscErrorCode ierr;

    if (J >= 0 && J < N &&
        K >= 0 && K < N) {
        ierr = MatSetValues (aimplfd->A, 1, &J, 1, &K, &C, INSERT_VALUES); CHKERR;
    }
}

sf_petsc_aimplfd2 sf_petsc_aimplfd2_init (int nz, int nx, float dz, float dx, float dt,
                                          float *v, int niter)
/*< Initialize implicit F-D object. >*/
{
    PetscInt J, K, N;
    int ix, iz, idx;
    float v2, rtx2, rtz2;
    PetscErrorCode ierr;
    PetscScalar Rox, Roz, Val;

    sf_petsc_aimplfd2 aimplfd;

    aimplfd = (sf_petsc_aimplfd2)sf_alloc (1, sizeof (struct PETScAimplFD2));

    aimplfd->Nx = nx;
    aimplfd->Nz = nz;
    aimplfd->dz = dz;
    aimplfd->dx = dx;
    aimplfd->dt = dt;
    aimplfd->comm = MPI_COMM_WORLD;

    /* Matrix */
    ierr = MatCreate (aimplfd->comm, &aimplfd->A); CHKERR;
    /* Matrix width/heigh */
    N = aimplfd->Nx*aimplfd->Nz;
    ierr = MatSetSizes (aimplfd->A, PETSC_DECIDE, PETSC_DECIDE,
                        N, N); CHKERR;
    ierr = MatSetType (aimplfd->A, MATMPIAIJ); CHKERR;
    /* Solver type */
    ierr = KSPCreate (aimplfd->comm, &aimplfd->Solver); CHKERR;
    ierr = KSPSetOperators (aimplfd->Solver, aimplfd->A, aimplfd->A, SAME_NONZERO_PATTERN); CHKERR;
    ierr = KSPSetType (aimplfd->Solver, KSPGMRES); CHKERR;
    ierr = KSPSetTolerances (aimplfd->Solver, PETSC_DEFAULT, PETSC_DEFAULT,
                             PETSC_DEFAULT, niter); CHKERR;
    /* Solution and r.h.s. vectors */
    ierr = VecCreateMPI (aimplfd->comm, PETSC_DECIDE, N, &aimplfd->Rhs); CHKERR;
    Val = 0.0;
    ierr = VecSet (aimplfd->Rhs, Val); CHKERR;
    ierr = VecDuplicate (aimplfd->Rhs, &aimplfd->Ut); CHKERR;
    ierr = VecDuplicate (aimplfd->Rhs, &aimplfd->Ut1); CHKERR;
    ierr = VecDuplicate (aimplfd->Rhs, &aimplfd->Ut2); CHKERR;

    /* Build matrix */
    ierr = MatAssemblyBegin (aimplfd->A, MAT_FINAL_ASSEMBLY); CHKERR;
    rtx2 = (dt*dt)/(dx*dx);
    rtz2 = (dt*dt)/(dz*dz);
    for (ix = 0; ix < nx; ix++) {
        for (iz = 0; iz < nz; iz++) {
            idx = ix*nz + iz;
            v2 = v[idx];
            v2 *= v2;
            Rox = v2*rtx2;
            Roz = v2*rtz2;
            J = idx;
            K = J - nz;
            sf_petsc_aimplfd2_set_mat (aimplfd, J, K, N, -Rox);
            K = J - 1;
            sf_petsc_aimplfd2_set_mat (aimplfd, J, K, N, -Roz);
            K = J;
            sf_petsc_aimplfd2_set_mat (aimplfd, J, J, N, 1 + 2*(Rox + Roz));
            K = J + 1;
            sf_petsc_aimplfd2_set_mat (aimplfd, J, K, N, -Roz);
            K = J + nz;
            sf_petsc_aimplfd2_set_mat (aimplfd, J, K, N, -Rox);
        }
    }
    /* Finish matrix assembly */
    ierr = MatAssemblyEnd (aimplfd->A, MAT_FINAL_ASSEMBLY); CHKERR;
    ierr = MatSetOption (aimplfd->A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE); CHKERR;

    return aimplfd;
}

void sf_petsc_aimplfd2_next_step (sf_petsc_aimplfd2 aimplfd)
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

void sf_petsc_aimplfd2_add_source (sf_petsc_aimplfd2 aimplfd, float f, int iz, int ix)
/*< Inject source at a specific location. >*/
{
    PetscInt J;
    PetscErrorCode ierr;
    PetscScalar Val;

    ierr = VecAssemblyBegin (aimplfd->Ut1); CHKERR;
    Val = f*0.5*aimplfd->dt*aimplfd->dt;
    J = aimplfd->Nz*ix + iz;
    ierr = VecSetValue (aimplfd->Ut1, J, Val, ADD_VALUES); CHKERR;
    ierr = VecAssemblyEnd (aimplfd->Ut1); CHKERR;
}

void sf_petsc_aimplfd2_get_wavefield (sf_petsc_aimplfd2 aimplfd, float *u)
/*< Get wavefield values at the current time step. >*/
{
    int i, nx, nz, ix, iz;
    PetscInt J;
    PetscErrorCode ierr;
    PetscScalar Val;
    Vec Uout;
    VecScatter Uctx;

    ierr = VecScatterCreateToZero (aimplfd->Ut, &Uctx, &Uout); CHKERR;
    ierr = VecScatterBegin (Uctx, aimplfd->Ut, Uout, INSERT_VALUES, SCATTER_FORWARD); CHKERR;
    nx = aimplfd->Nx;
    nz = aimplfd->Nz;
    for (ix = 0; ix < nx; ix++) {
        for (iz = 0; iz < nz; iz++) {
            J = ix*nz + iz;
            i = J;
            ierr = VecGetValues (Uout, 1, &J, &Val); CHKERR;
            u[i] = Val;
        }
    }
    ierr = VecScatterEnd (Uctx, aimplfd->Ut, Uout, INSERT_VALUES, SCATTER_FORWARD); CHKERR;
    ierr = VecScatterDestroy (Uctx); CHKERR;
    ierr = VecDestroy (Uout); CHKERR;
}

void sf_petsc_aimplfd2_destroy (sf_petsc_aimplfd2 aimplfd)
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
