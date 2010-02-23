/* Implicit solution of 1-D acoustic wave equation. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
#include <petscksp.h>

int main (int argc, char* argv[]) {
    int nx, nt, i, it, is;
    float ox, ot, dx, dt;
    float ro, sx, rtx2;
    float *u, *v, *f;
    /* I/O */
    sf_file usol, vel, src;
    /* PETSc */
    PetscInt J, K, N;
    MPI_Comm comm;
    PetscErrorCode ierr;
    KSP Solver;
    Mat A;
    Vec Rhs, Ut, Ut1, Ut2, Uout;
    VecScatter Uctx;
    PetscScalar Ro12, Rom, Val;
    PetscInt Niter;

    sf_init (argc, argv);

    vel = sf_input ("in");
    /* Velocity */
    usol = sf_output ("out");
    /* Solution - wavefield */

    if (SF_FLOAT != sf_gettype (vel))
        sf_error ("Need float input");
    if (!sf_histint (vel, "n1", &nx)) sf_error ("No n1= in input");
    if (!sf_histfloat (vel, "d1", &dx)) sf_error ("No d1= in input");
    if (!sf_histfloat (vel, "o1", &ox)) ox = 0.;

    if (!sf_getstring ("src")) sf_error ("Need src=");
    /* Source wavelet */
    src = sf_input ("src");
    if (!sf_histint (src, "n1", &nt)) sf_error ("No n1= in src");
    if (!sf_histfloat (src, "d1", &dt)) sf_error ("No d1= in src");
    if (!sf_histfloat (src, "o1", &ot)) ot = 0.;

    if (!sf_histfloat (src, "o2", &sx)) sf_error ("No o2= in src");
    is = (sx - ox)/dx;
    if (is < 0 || is >= nx) sf_error ("Shot location is out of bounds");

    if (!sf_getint ("niter", &i)) i = 10;
    /* Number of solver iterations */
    Niter = i;

    /* Wavefield */
    u = sf_floatalloc (nx);

    /* Velocity */
    v = sf_floatalloc (nx);
    /* Source */
    f = sf_floatalloc (nt);

    /* Read velocity */
    sf_floatread (v, nx, vel);
    /* Read source */
    sf_floatread (f, nt, src);

    /* Create time axis in output */
    sf_putint (usol, "n2", nt);
    sf_putfloat (usol, "d2", dt);
    sf_putfloat (usol, "o2", ot);

    /********************************
     *         PETSc part           *
     ********************************/

    /* Initialization */
    ierr = PetscInitialize (&argc, &argv, 0, 0); CHKERRQ(ierr);
    comm = MPI_COMM_WORLD;

    PetscFPrintf (comm, stderr, "\nCreating PETSc data structures\n");

    /* Matrix */
    ierr = MatCreate (comm, &A); CHKERRQ(ierr);
    /* Matrix width/heigh */
    N = nx;
    ierr = MatSetSizes (A, PETSC_DECIDE, PETSC_DECIDE,
                        N, N); CHKERRQ(ierr);
    ierr = MatSetType (A, MATMPIAIJ); CHKERRQ(ierr);
    /* Solver type */
    ierr = KSPCreate (comm, &Solver); CHKERRQ(ierr);
    ierr = KSPSetOperators (Solver, A, A, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
    ierr = KSPSetType (Solver, KSPGMRES); CHKERRQ(ierr);
    ierr = KSPSetTolerances (Solver, PETSC_DEFAULT, PETSC_DEFAULT,
                             PETSC_DEFAULT, Niter); CHKERRQ(ierr);
    /* Solution and r.h.s. vectors */
    ierr = VecCreateMPI (comm, PETSC_DECIDE, N, &Rhs); CHKERRQ(ierr);
    ierr = VecDuplicate (Rhs, &Ut); CHKERRQ(ierr);
    ierr = VecDuplicate (Rhs, &Ut1); CHKERRQ(ierr);
    ierr = VecDuplicate (Rhs, &Ut2); CHKERRQ(ierr);

    PetscFPrintf (comm, stderr, "Preparing PETSc data\n");

    /* Build matrix */
    ierr = MatAssemblyBegin (A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    rtx2 = (dt*dt)/(dx*dx);
    /* Left B.C. */
    J = 0;
    K = 0;
    ro = v[0]*rtx2;
    Ro12 = 1 + 2*ro;
    Rom = -ro;
    ierr = MatSetValues (A, 1, &J, 1, &K, &Ro12, INSERT_VALUES); CHKERRQ(ierr);
    K = K + 1;
    ierr = MatSetValues (A, 1, &J, 1, &K, &Rom, INSERT_VALUES); CHKERRQ(ierr);
    for (i = 1; i < (nx - 1); i++) {
        J = i;
        K = i - 1;
        ro = v[i]*rtx2;
        /* Coefficients */
        Ro12 = 1 + 2*ro;
        Rom = -ro;
        ierr = MatSetValues (A, 1, &J, 1, &K, &Rom, INSERT_VALUES); CHKERRQ(ierr);
        K = K + 1;
        ierr = MatSetValues (A, 1, &J, 1, &K, &Ro12, INSERT_VALUES); CHKERRQ(ierr);
        K = K + 1;
        ierr = MatSetValues (A, 1, &J, 1, &K, &Rom, INSERT_VALUES); CHKERRQ(ierr);
    }
    /* Right B.C. */
    J = nx - 1;
    K = nx - 2;
    ro = v[nx - 1]*rtx2;
    Ro12 = 1 + 2*ro;
    Rom = -ro;
    ierr = MatSetValues (A, 1, &J, 1, &K, &Rom, INSERT_VALUES); CHKERRQ(ierr);
    K = K + 1;
    ierr = MatSetValues (A, 1, &J, 1, &K, &Ro12, INSERT_VALUES); CHKERRQ(ierr);
    /* Finish matrix assembly */
    ierr = MatAssemblyEnd (A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatSetOption (A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE); CHKERRQ(ierr);

    PetscFPrintf (comm, stderr, "Running GMRES solver\n");

    /* Loop in time */
    for (it = 0; it < nt; it++) {
        PetscFPrintf (comm, stderr, "Timestep #%d, t=%f\n", it, it*dt);

        /* R.H.S. */
        ierr = VecAssemblyBegin (Rhs); CHKERRQ(ierr);
        Val = -1.0;
        VecScale (Ut2, Val);
        Val = 2.0;
        VecWAXPY (Rhs, Val, Ut1, Ut2);
        ierr = VecAssemblyEnd (Rhs); CHKERRQ(ierr);

        /* Solve the system */
        ierr = KSPSolve (Solver, Rhs, Ut); CHKERRQ(ierr);

        /* Move time steps */
        ierr = VecCopy (Ut1, Ut2); CHKERRQ(ierr);
        ierr = VecCopy (Ut, Ut1); CHKERRQ(ierr);

        /* Insert source */
        ierr = VecAssemblyBegin (Ut1); CHKERRQ(ierr);
        Val = f[it]*0.5*dt*dt;
        J = is;
        ierr = VecSetValue (Ut1, J, Val, ADD_VALUES); CHKERRQ(ierr);
        ierr = VecAssemblyEnd (Ut1); CHKERRQ(ierr);

        /* Get solution from PETSc */
        VecScatterCreateToZero (Ut, &Uctx, &Uout);
        VecScatterBegin (Uctx, Ut, Uout, INSERT_VALUES, SCATTER_FORWARD);
        for (i = 0; i < nx; i++) {
            J = i;
            VecGetValues (Uout, 1, &J, &Val);
            u[i] = Val;
        }
        VecScatterEnd (Uctx, Ut, Uout, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterDestroy (Uctx);
        VecDestroy (Uout);

        /* Write the solution */
        sf_floatwrite (u, nx, usol);
    }

    /* Clean up */
    ierr = KSPDestroy (Solver); CHKERRQ(ierr);
    ierr = VecDestroy (Rhs); CHKERRQ(ierr);
    ierr = VecDestroy (Ut); CHKERRQ(ierr);
    ierr = VecDestroy (Ut1); CHKERRQ(ierr);
    ierr = VecDestroy (Ut2); CHKERRQ(ierr);
    ierr = PetscFinalize ();

    exit (0);
}
