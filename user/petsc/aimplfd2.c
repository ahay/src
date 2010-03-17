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
    PetscInt Nx, Nz, Npad, Nxpad, Nzpad;
    MPI_Comm comm;
    KSP Solver;
    Mat A;
    Vec Ut, Ut1, Ut2;
    float dz, dx, dt;
};
/* concrete data type */

#define SCHKERR CHKERRABORT(comm,ierr)
#define CHKERR CHKERRABORT(aimplfd->comm,ierr)

/* Show contents of a matrix */
static void sf_petsc_mat_view (MPI_Comm comm, Mat A) {
    PetscViewer viewer;
    PetscErrorCode ierr;

/*
    ierr = PetscViewerDrawOpen (comm, NULL, "Matrix", PETSC_DECIDE, PETSC_DECIDE,
                                                      PETSC_DECIDE, PETSC_DECIDE, &viewer); SCHKERR;
*/
    ierr = PetscViewerASCIIOpen (comm, "mat.out", &viewer); SCHKERR;
    ierr = MatView (A, viewer); SCHKERR;
    ierr = PetscViewerDestroy (viewer); SCHKERR;
}
    
/* Remove padding from indices and get velocity from the array */
static float sf_petsc_aimplfd2_get_vel (sf_petsc_aimplfd2 aimplfd, float *v, int ix, int iz) {
    int ixs, izs;
    ixs = ix - aimplfd->Npad;
    izs = iz - aimplfd->Npad;
    if (ix < aimplfd->Npad)
        ixs = 0;
    else if (ix >= (aimplfd->Npad + aimplfd->Nx))
        ixs = aimplfd->Nx - 1;
    if (iz < aimplfd->Npad)
        izs = 0;
    else if (iz >= (aimplfd->Npad + aimplfd->Nz))
        izs = aimplfd->Nz - 1;
    return v[ixs*aimplfd->Nz + izs];
}

/* Check boundaries and insert a coefficient into the matrix */
static void sf_petsc_aimplfd2_set_mat (MPI_Comm comm, Mat A, PetscInt J, PetscInt K, PetscInt N, PetscScalar C)
{
    PetscErrorCode ierr;

    if (J >= 0 && J < N &&
        K >= 0 && K < N) {
        ierr = MatSetValues (A, 1, &J, 1, &K, &C, INSERT_VALUES); SCHKERR;
    }
}

sf_petsc_aimplfd2 sf_petsc_aimplfd2_init (int nz, int nx, float dz, float dx, float dt,
                                          float *v, int niter)
/*< Initialize implicit F-D object. >*/
{
    PetscInt J, K, N;
    int ix, iz;
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
    aimplfd->Npad = 2; /* Pad area thickness */
    aimplfd->Nxpad = aimplfd->Nx + aimplfd->Npad*2; /* Total thickness in x */
    aimplfd->Nzpad = aimplfd->Nz + aimplfd->Npad*2; /* Total thickness in x */
    aimplfd->comm = MPI_COMM_WORLD;

    /* Matrix */
    ierr = MatCreate (aimplfd->comm, &aimplfd->A); CHKERR;
    /* Matrix width/heigh */
    N = aimplfd->Nxpad*aimplfd->Nzpad;
    ierr = MatSetSizes (aimplfd->A, PETSC_DECIDE, PETSC_DECIDE,
                        N, N); CHKERR;
    ierr = MatSetType (aimplfd->A, MATMPIAIJ); CHKERR;
    ierr = MatSetOption (aimplfd->A, MAT_SYMMETRIC, PETSC_TRUE); CHKERR;
    /* Solver type */
    ierr = KSPCreate (aimplfd->comm, &aimplfd->Solver); CHKERR;
    ierr = KSPSetOperators (aimplfd->Solver, aimplfd->A, aimplfd->A, SAME_NONZERO_PATTERN); CHKERR;
    ierr = KSPSetType (aimplfd->Solver, KSPCGS); CHKERR;
    ierr = KSPSetTolerances (aimplfd->Solver, PETSC_DEFAULT, PETSC_DEFAULT,
                             PETSC_DEFAULT, niter); CHKERR;
    /* Solution and r.h.s. vectors */
    ierr = VecCreateMPI (aimplfd->comm, PETSC_DECIDE, N, &aimplfd->Ut); CHKERR;
    ierr = VecZeroEntries (aimplfd->Ut); CHKERR;
    ierr = VecDuplicate (aimplfd->Ut, &aimplfd->Ut1); CHKERR;
    ierr = VecDuplicate (aimplfd->Ut, &aimplfd->Ut2); CHKERR;

    /* Build propagation matrix */
    ierr = MatAssemblyBegin (aimplfd->A, MAT_FINAL_ASSEMBLY); CHKERR;
    rtx2 = (dt*dt)/(dx*dx);
    rtz2 = (dt*dt)/(dz*dz);
    nx = aimplfd->Nxpad;
    nz = aimplfd->Nzpad;
    Val = 1.0;
    for (ix = 0; ix < nx; ix++) {
        for (iz = 0; iz < nz; iz++) {
            if (ix < aimplfd->Npad || ix >= (aimplfd->Nx + aimplfd->Npad) ||
                iz < aimplfd->Npad || iz >= (aimplfd->Nz + aimplfd->Npad)) {
                J = ix*nz + iz;
                sf_petsc_aimplfd2_set_mat (aimplfd->comm, aimplfd->A, J, J, N, Val);
            } else {
                v2 = sf_petsc_aimplfd2_get_vel (aimplfd, v, ix, iz);
                v2 *= v2;
                Rox = v2*rtx2;
                Roz = v2*rtz2;
                J = ix*nz + iz;
                K = J - nz;
                sf_petsc_aimplfd2_set_mat (aimplfd->comm, aimplfd->A, J, K, N, -Rox);
                K = J - 1;
                sf_petsc_aimplfd2_set_mat (aimplfd->comm, aimplfd->A, J, K, N, -Roz);
                K = J;
                sf_petsc_aimplfd2_set_mat (aimplfd->comm, aimplfd->A, J, J, N, 1 + 2*(Rox + Roz));
                K = J + 1;
                sf_petsc_aimplfd2_set_mat (aimplfd->comm, aimplfd->A, J, K, N, -Roz);
                K = J + nz;
                sf_petsc_aimplfd2_set_mat (aimplfd->comm, aimplfd->A, J, K, N, -Rox);
            }
        }
    }
    /* Finish matrix assembly */
    ierr = MatAssemblyEnd (aimplfd->A, MAT_FINAL_ASSEMBLY); CHKERR;
    ierr = MatSetOption (aimplfd->A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE); CHKERR;
/*
    sf_petsc_mat_view (aimplfd->comm, aimplfd->A);
*/
    return aimplfd;
}

void sf_petsc_aimplfd2_next_step (sf_petsc_aimplfd2 aimplfd)
/*< Calculate next time step. >*/
{
    PetscErrorCode ierr;
    PetscScalar Val1, Val2;
    PetscInt N;

    /* R.H.S. */
    Val1 = 2.0;
    Val2 = -1.0;
    /* Ut2 becomes right-hand side */
    VecAXPBY (aimplfd->Ut2, Val1, Val2, aimplfd->Ut1);

    /* Solve the system */
    ierr = KSPSolve (aimplfd->Solver, aimplfd->Ut2, aimplfd->Ut); CHKERR;
    ierr = KSPGetIterationNumber (aimplfd->Solver, &N); CHKERR;
    PetscFPrintf (MPI_COMM_WORLD, stderr, "Iterations number - %d\n", N);
    /* Move time steps */
    ierr = VecCopy (aimplfd->Ut1, aimplfd->Ut2); CHKERR;
    ierr = VecCopy (aimplfd->Ut, aimplfd->Ut1); CHKERR;
}

int sf_petsc_aimplfd2_get_niter (sf_petsc_aimplfd2 aimplfd)
/*< Get number of iterations from the previous solver run. >*/
{
    PetscErrorCode ierr;
    PetscInt N;

    ierr = KSPGetIterationNumber (aimplfd->Solver, &N); CHKERR;

    return N;
}

void sf_petsc_aimplfd2_add_source_ut1 (sf_petsc_aimplfd2 aimplfd, float f, int iz, int ix)
/*< Inject source at a specific location into t-1 step. >*/
{
    PetscInt J;
    PetscErrorCode ierr;
    PetscScalar Val;

    ierr = VecAssemblyBegin (aimplfd->Ut1); CHKERR;
    Val = f*0.5*aimplfd->dt*aimplfd->dt;
    J = aimplfd->Nzpad*(ix + aimplfd->Npad) + iz + aimplfd->Npad;
    ierr = VecSetValue (aimplfd->Ut1, J, Val, ADD_VALUES); CHKERR;
    ierr = VecAssemblyEnd (aimplfd->Ut1); CHKERR;
}

void sf_petsc_aimplfd2_add_source_ut2 (sf_petsc_aimplfd2 aimplfd, float f, int iz, int ix)
/*< Inject source at a specific location into t-1 step. >*/
{
    PetscInt J;
    PetscErrorCode ierr;
    PetscScalar Val;

    ierr = VecAssemblyBegin (aimplfd->Ut2); CHKERR;
    Val = f*0.5*aimplfd->dt*aimplfd->dt;
    J = aimplfd->Nzpad*(ix + aimplfd->Npad) + iz + aimplfd->Npad;
    ierr = VecSetValue (aimplfd->Ut2, J, Val, ADD_VALUES); CHKERR;
    ierr = VecAssemblyEnd (aimplfd->Ut2); CHKERR;
}

void sf_petsc_aimplfd2_get_wavefield_ut1 (sf_petsc_aimplfd2 aimplfd, float *u)
/*< Get wavefield values at the current time step. >*/
{
    int i, nx, nz, ix, iz;
    PetscInt J;
    PetscErrorCode ierr;
    PetscScalar Val;
    Vec Uout;
    VecScatter Uctx;

    ierr = VecScatterCreateToZero (aimplfd->Ut1, &Uctx, &Uout); CHKERR;
    ierr = VecScatterBegin (Uctx, aimplfd->Ut1, Uout, INSERT_VALUES, SCATTER_FORWARD); CHKERR;
    nx = aimplfd->Nx + aimplfd->Npad;
    nz = aimplfd->Nz + aimplfd->Npad;
    i = 0;
    for (ix = aimplfd->Npad; ix < nx; ix++) {
        for (iz = aimplfd->Npad; iz < nz; iz++) {
            J = ix*aimplfd->Nzpad + iz;
            ierr = VecGetValues (Uout, 1, &J, &Val); CHKERR;
            u[i] = Val;
            i++;
        }
    }
    ierr = VecScatterEnd (Uctx, aimplfd->Ut1, Uout, INSERT_VALUES, SCATTER_FORWARD); CHKERR;
    ierr = VecScatterDestroy (Uctx); CHKERR;
    ierr = VecDestroy (Uout); CHKERR;
}

void sf_petsc_aimplfd2_get_wavefield_ut2 (sf_petsc_aimplfd2 aimplfd, float *u)
/*< Get wavefield values at the current time step. >*/
{
    int i, nx, nz, ix, iz;
    PetscInt J;
    PetscErrorCode ierr;
    PetscScalar Val;
    Vec Uout;
    VecScatter Uctx;

    ierr = VecScatterCreateToZero (aimplfd->Ut2, &Uctx, &Uout); CHKERR;
    ierr = VecScatterBegin (Uctx, aimplfd->Ut2, Uout, INSERT_VALUES, SCATTER_FORWARD); CHKERR;
    nx = aimplfd->Nx + aimplfd->Npad;
    nz = aimplfd->Nz + aimplfd->Npad;
    i = 0;
    for (ix = aimplfd->Npad; ix < nx; ix++) {
        for (iz = aimplfd->Npad; iz < nz; iz++) {
            J = ix*aimplfd->Nzpad + iz;
            ierr = VecGetValues (Uout, 1, &J, &Val); CHKERR;
            u[i] = Val;
            i++;
        }
    }
    ierr = VecScatterEnd (Uctx, aimplfd->Ut2, Uout, INSERT_VALUES, SCATTER_FORWARD); CHKERR;
    ierr = VecScatterDestroy (Uctx); CHKERR;
    ierr = VecDestroy (Uout); CHKERR;
}

void sf_petsc_aimplfd2_destroy (sf_petsc_aimplfd2 aimplfd)
/*< Destroy the object. >*/
{
    PetscErrorCode ierr;

    ierr = KSPDestroy (aimplfd->Solver); CHKERR;
    ierr = MatDestroy (aimplfd->A); CHKERR;
    ierr = VecDestroy (aimplfd->Ut); CHKERR;
    ierr = VecDestroy (aimplfd->Ut1); CHKERR;
    ierr = VecDestroy (aimplfd->Ut2); CHKERR;   

    free (aimplfd);
}

