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

#include <assert.h>

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
    MatNullSpace Nullsp;
    Vec Ut, Ut1, Ut2;
    float dz, dx, dt;
    float *bzl, *bzh, *bxl, *bxh; /* B.C. */
};
/* concrete data type */

#define SCHKERR CHKERRABORT(comm,ierr)
#define CHKERR CHKERRABORT(aimplfd->comm,ierr)

/* Show contents of a matrix */
void sf_petsc_mat_view (MPI_Comm comm, Mat A) {
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

#define MAX_RC 9

/* Null space remover */
static PetscErrorCode sf_petsc_aimplfd2_null_space (Vec x, void *ctx) {
    sf_petsc_aimplfd2 aimplfd = (sf_petsc_aimplfd2)ctx;
    PetscErrorCode ierr;
    PetscInt idx;
    PetscScalar Val = 0.0;
    int nx, nz, ix, iz;

    nx = aimplfd->Nxpad;
    nz = aimplfd->Nzpad;
    /* Cleanup the padded area */
    ierr = VecAssemblyBegin (x); CHKERR;
    for (ix = 0; ix < nx; ix++) {
        for (iz = 0; iz < nz; iz++) {
            idx = ix*nz + iz;
            if (ix < aimplfd->Npad || ix >= (aimplfd->Npad + aimplfd->Nx) ||
                iz < aimplfd->Npad || iz >= (aimplfd->Npad + aimplfd->Nz)) {
                ierr = VecSetValue (x, idx, Val, INSERT_VALUES); CHKERR;
            }
        }
    }
    ierr = VecAssemblyEnd (x); CHKERR;

    PetscFunctionReturn (0);
}

sf_petsc_aimplfd2 sf_petsc_aimplfd2_init (int nz, int nx, float dz, float dx, float dt,
                                          float *v, int niter, bool fourth)
/*< Initialize implicit F-D object. >*/
{
    PetscInt J, N;
    int ix, iz, idx;
    float v2, rtx2, rtz2, d;
    PetscErrorCode ierr;
    PetscScalar Rox, Roz, Val, Vals[MAX_RC];
    PetscInt Cols[MAX_RC];

    sf_petsc_aimplfd2 aimplfd;

    aimplfd = (sf_petsc_aimplfd2)sf_alloc (1, sizeof (struct PETScAimplFD2));

    aimplfd->Nx = nx;
    aimplfd->Nz = nz;
    aimplfd->dz = dz;
    aimplfd->dx = dx;
    aimplfd->dt = dt;
    aimplfd->Npad = 5; /* Pad area thickness */
    aimplfd->Nxpad = aimplfd->Nx + aimplfd->Npad*2; 
    aimplfd->Nzpad = aimplfd->Nz + aimplfd->Npad*2; 
    aimplfd->comm = MPI_COMM_WORLD;

    /* Matrix width/heigh */
    N = aimplfd->Nxpad*aimplfd->Nzpad;

    /* Matrix */
    /*
    ierr = MatCreate (aimplfd->comm, &aimplfd->A); CHKERR;
    ierr = MatSetSizes (aimplfd->A, PETSC_DECIDE, PETSC_DECIDE,
                        N, N); CHKERR;
    ierr = MatSetType (aimplfd->A, MATMPIBAIJ); CHKERR;
    */
    ierr = MatCreateSeqAIJ (aimplfd->comm, N, N, fourth ? 9 : 5, PETSC_NULL, &aimplfd->A); CHKERR;
    ierr = MatZeroEntries (aimplfd->A); CHKERR;
    ierr = MatSetOption (aimplfd->A, MAT_SYMMETRIC, PETSC_TRUE); CHKERR;
    /*
    ierr = MatSetOption (aimplfd->A, MAT_USE_HASH_TABLE, PETSC_TRUE); CHKERR;
    */
    /* Solution and r.h.s. vectors */
    ierr = VecCreateMPI (aimplfd->comm, PETSC_DECIDE, N, &aimplfd->Ut); CHKERR;
    ierr = VecZeroEntries (aimplfd->Ut); CHKERR;
    ierr = VecDuplicate (aimplfd->Ut, &aimplfd->Ut1); CHKERR;
    ierr = VecDuplicate (aimplfd->Ut, &aimplfd->Ut2); CHKERR;

    /* Build propagation matrix */
    rtx2 = (dt*dt)/(dx*dx);
    rtz2 = (dt*dt)/(dz*dz);
    nx = aimplfd->Nxpad;
    nz = aimplfd->Nzpad;
    Val = 1.0;
    for (ix = 0; ix < nx; ix++) {
        for (iz = 0; iz < nz; iz++) {
            idx = ix*nz + iz;
            if (ix < aimplfd->Npad || ix >= (aimplfd->Npad + aimplfd->Nx) ||
                iz < aimplfd->Npad || iz >= (aimplfd->Npad + aimplfd->Nz)) {
                J = idx;
                ierr = MatSetValue (aimplfd->A, J, J, Val, INSERT_VALUES); CHKERR;
            } else {
                v2 = sf_petsc_aimplfd2_get_vel (aimplfd, v, ix, iz);
                v2 *= v2;
                Rox = v2*rtx2;
                Roz = v2*rtz2;
                if (fourth) {
                    J = idx;
                    Cols[0] = J - 2*nz;
                    Cols[1] = J - nz;
                    Cols[2] = J - 2;
                    Cols[3] = J - 1;
                    Cols[4] = J;
                    Cols[5] = J + 1;
                    Cols[6] = J + 2;
                    Cols[7] = J + nz;
                    Cols[8] = J + 2*nz;
                    Vals[0] = 1.0/12.0*Rox;
                    Vals[1] = -16.0/12.0*Rox;
                    Vals[2] = 1.0/12.0*Roz;
                    Vals[3] = -16.0/12.0*Roz;
                    Vals[4] = 1 + 30.0/12.0*(Rox + Roz);
                    Vals[5] = -16.0/12.0*Roz;
                    Vals[6] = 1.0/12.0*Roz;
                    Vals[7] = -16.0/12.0*Rox;
                    Vals[8] = 1.0/12.0*Rox;
                    ierr = MatSetValues (aimplfd->A, 1, &J, 9, (const PetscInt *)&Cols, Vals, INSERT_VALUES); CHKERR;
                } else {
                    J = idx;
                    Cols[0] = J - nz;
                    Cols[1] = J - 1;
                    Cols[2] = J;
                    Cols[3] = J + 1;
                    Cols[4] = J + nz;
                    Vals[0] = -Rox;
                    Vals[1] = -Roz;
                    Vals[2] = 1 + 2.0*(Rox + Roz);
                    Vals[3] = -Roz;
                    Vals[4] = -Rox;
                    ierr = MatSetValues (aimplfd->A, 1, &J, 5,  (const PetscInt *)&Cols, Vals, INSERT_VALUES); CHKERR;
                }
            }
        }
    }
    /* Finish matrix assembly */
    ierr = MatAssemblyBegin (aimplfd->A, MAT_FINAL_ASSEMBLY); CHKERR;
    ierr = MatAssemblyEnd (aimplfd->A, MAT_FINAL_ASSEMBLY); CHKERR;
    ierr = MatSetOption (aimplfd->A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE); CHKERR;
/*
    sf_petsc_mat_view (aimplfd->comm, aimplfd->A);
*/

    aimplfd->bzl = sf_floatalloc (aimplfd->Nxpad);
    aimplfd->bzh = sf_floatalloc (aimplfd->Nxpad);
    aimplfd->bxl = sf_floatalloc (aimplfd->Nzpad);
    aimplfd->bxh = sf_floatalloc (aimplfd->Nzpad);

    rtx2 = dt/dx;
    rtz2 = dt/dz;
    /* Top/Bottom coefficients for B.C. */
    for (ix = 0; ix < aimplfd->Nxpad; ix++) {
        v2 = sf_petsc_aimplfd2_get_vel (aimplfd, v, ix, 0);
        d = v2*rtz2;
        aimplfd->bzl[ix] = (1-d)/(1+d);
        v2 = sf_petsc_aimplfd2_get_vel (aimplfd, v, ix, aimplfd->Nzpad);
        d = v2*rtz2;
        aimplfd->bzh[ix] = (1-d)/(1+d);
    }
    /* Left/right coefficients for B.C. */
    for (iz = 0; iz < aimplfd->Nzpad; iz++) {
        v2 = sf_petsc_aimplfd2_get_vel (aimplfd, v, 0, iz);
        d = v2*rtx2;
        aimplfd->bxl[iz] = (1-d)/(1+d);
        v2 = sf_petsc_aimplfd2_get_vel (aimplfd, v, aimplfd->Nxpad, iz);
        d = v2*rtx2;
        aimplfd->bxh[iz] = (1-d)/(1+d);
    }

    ierr = MatNullSpaceCreate (aimplfd->comm, PETSC_FALSE, 0, PETSC_NULL, &aimplfd->Nullsp); CHKERR;
    ierr = MatNullSpaceSetFunction (aimplfd->Nullsp, sf_petsc_aimplfd2_null_space, aimplfd); CHKERR;

    /* Solver type */
    ierr = KSPCreate (aimplfd->comm, &aimplfd->Solver); CHKERR;
    ierr = KSPSetOperators (aimplfd->Solver, aimplfd->A, aimplfd->A, DIFFERENT_NONZERO_PATTERN); CHKERR;
    ierr = KSPSetType (aimplfd->Solver, KSPGMRES); CHKERR;
    ierr = KSPSetTolerances (aimplfd->Solver, PETSC_DEFAULT, PETSC_DEFAULT,
                             PETSC_DEFAULT, niter); CHKERR;
    ierr = KSPSetNullSpace (aimplfd->Solver, aimplfd->Nullsp); CHKERR;
    ierr = KSPSetUp (aimplfd->Solver); CHKERR;

    return aimplfd;
}

/* Apply B.C. */
static void sf_petsc_aimplfd2_bc_apply (sf_petsc_aimplfd2 aimplfd) {
    int iz, ix, iop;
    int nop = aimplfd->Npad;
    PetscErrorCode ierr;
    PetscScalar **uo, **um;

    ierr = VecGetArray2d (aimplfd->Ut1, aimplfd->Nxpad, aimplfd->Nzpad, 0, 0, &um); CHKERR;
    ierr = VecGetArray2d (aimplfd->Ut2, aimplfd->Nxpad, aimplfd->Nzpad, 0, 0, &uo); CHKERR;

    for (ix = 0; ix < aimplfd->Nxpad; ix++) {
        for (iop = 0; iop < nop; iop++) {
            /* top BC */
            iz = nop - iop;
            uo[ix][iz] = um[ix][iz + 1] +(um[ix][iz] - uo[ix][iz + 1])*aimplfd->bzl[ix];
            /* bottom BC */
            iz = aimplfd->Nzpad - nop + iop - 1;
            uo[ix][iz] = um[ix][iz - 1] +(um[ix][iz] - uo[ix][iz - 1])*aimplfd->bzh[ix];
        }
    }

    for (iz = 0; iz < aimplfd->Nzpad; iz++) {
        for (iop = 0; iop < nop; iop++) {
            /* left BC */
            ix = nop - iop;
            uo[ix][iz] = um[ix + 1][iz] + (um[ix][iz] - uo[ix + 1][iz])*aimplfd->bxl[iz];
            /* right BC */
            ix = aimplfd->Nxpad - nop + iop - 1;
            uo[ix  ][iz] = um[ix - 1][iz] + (um[ix][iz] - uo[ix - 1][iz])*aimplfd->bxh[iz];
        }
    }

    ierr = VecRestoreArray2d (aimplfd->Ut2, aimplfd->Nzpad, aimplfd->Nxpad, 0, 0, &uo); CHKERR;
    ierr = VecRestoreArray2d (aimplfd->Ut1, aimplfd->Nzpad, aimplfd->Nxpad, 0, 0, &um); CHKERR;
    return;
}

void sf_petsc_aimplfd2_next_step (sf_petsc_aimplfd2 aimplfd)
/*< Calculate next time step. >*/
{
    PetscErrorCode ierr;
    PetscScalar Val1, Val2;

    /* R.H.S. */
    Val1 = 2.0;
    Val2 = -1.0;
    /* Ut2 becomes right-hand side */
    VecAXPBY (aimplfd->Ut2, Val1, Val2, aimplfd->Ut1);

    /* Apply B.C. */
    sf_petsc_aimplfd2_bc_apply (aimplfd);

    /* Solve the system */
    ierr = KSPSolve (aimplfd->Solver, aimplfd->Ut2, aimplfd->Ut); CHKERR;
    /*
    KSPConvergedReason reason;
    ierr = KSPGetConvergedReason (aimplfd->Solver, &reason); CHKERR;
    PetscFPrintf (MPI_COMM_WORLD, stderr, "Convergence reason: %d\n", reason);
    */

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

void sf_petsc_aimplfd2_add_source_ut2 (sf_petsc_aimplfd2 aimplfd, float f, int iz, int ix)
/*< Inject source at a specific location. >*/
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
    ierr = MatNullSpaceDestroy (aimplfd->Nullsp); CHKERR;

    free (aimplfd->bzl);
    free (aimplfd->bzh);
    free (aimplfd->bxl);
    free (aimplfd->bxh);

    free (aimplfd);
}

