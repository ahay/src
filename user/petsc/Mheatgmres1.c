/* Solution of 1-D heat equation with GMRES. */
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
    int nx, nt, i, j, k, it;
    float ox, dx, dt;
    float *T;
    float alpha, Ta, Tb, c;
    /* I/O */
    sf_file tinitial, solution;
    /* PETSc */
    PetscInt n;
    MPI_Comm comm;
    PetscErrorCode ierr;
    KSP Solver;
    Mat A;
    Vec Rhs, Sol;
    PetscScalar lc, cc, rc, ua, ub, t;
    PetscInt niter = 100;

    sf_init (argc, argv);

    tinitial = sf_input ("in");
    /* Initial temperature */
    solution = sf_output ("out");
    /* Solution - temperatures */

    if (SF_FLOAT != sf_gettype (tinitial))
        sf_error("Need float input");
    if (!sf_histint (tinitial, "n1", &nx)) sf_error ("No n1= in input");
    if (!sf_histfloat (tinitial, "d1", &dx)) sf_error ("No d1= in input");
    if (!sf_histfloat (tinitial, "o1", &ox)) ox = 0.;

    if (!sf_getfloat ("alpha", &alpha)) sf_error ("Need alpha=");
    /* Conductivity */

    if (!sf_getfloat ("dt", &dt)) sf_error ("Need dt=");
    /* Time step */
    if (!sf_getint ("nt", &nt)) sf_error ("Need nt=");
    /* Number of time steps */

    if (!sf_getfloat ("Ta", &Ta)) Ta = 0;
    /* Boundary condition on the left */
    if (!sf_getfloat ("Tb", &Tb)) Tb = 0;
    /* Boundary condition on the right */

    T  = sf_floatalloc (nx);
    /* Read the initial temperature */
    sf_floatread (T, nt, tinitial);

    /* Create time axis */
    sf_putint (solution, "n2", nt);
    sf_putfloat (solution, "o2", 0.);
    sf_putfloat (solution, "d2", dt);
    sf_putstring (solution, "label2", "Time");
    sf_putstring (solution, "unit2", "s");

    /********************************
     *         PETSc part           *
     ********************************/

    /* Initialization */
    ierr = PetscInitialize (&argc, &argv, 0, 0); CHKERRQ(ierr);
    comm = MPI_COMM_WORLD;
    /* Matrix */
    ierr = MatCreate (comm, &A); CHKERRQ(ierr);
    /* Matrix width/heigh */
    n = nx - 2; /* Two boundary points do not particpate */
    ierr = MatSetSizes (A, PETSC_DECIDE, PETSC_DECIDE,
                        n, n); CHKERRQ(ierr);
    ierr = MatSetType (A, MATMPIAIJ); CHKERRQ(ierr);
    /* Solver type */
    ierr = KSPCreate (comm, &Solver); CHKERRQ(ierr);
    ierr = KSPSetOperators (Solver, A, A, SAME_NONZERO_PATTERN); CHKERRQ(ierr);
    ierr = KSPSetType (Solver, KSPGMRES); CHKERRQ(ierr);
    ierr = KSPSetTolerances (Solver, PETSC_DEFAULT, PETSC_DEFAULT,
                             PETSC_DEFAULT, niter); CHKERRQ(ierr);
    /* Solution and r.h.s. vectors */
    ierr = VecCreateMPI (comm, PETSC_DECIDE, n, &Rhs); CHKERRQ(ierr);
    ierr = VecDuplicate (Rhs, &Sol); CHKERRQ(ierr);

    PetscFPrintf (comm, stderr, "\nRunning GMRES solver\n\n");

    /* Coefficients */
    c = dt*alpha / (dx*dx);
    lc = -c; /* Left of the diagonal */
    cc = 1 + 2.0*c; /* Diagonal */
    rc = -c; /* Right of the diagonal */
    /* Build matrix */
    ierr = MatAssemblyBegin (A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    for (i = 2; i < n; i++) {
        j = i - 1;
        k = i - 2;
        ierr = MatSetValues (A, 1, &j, 1, &k, &lc, INSERT_VALUES); CHKERRQ(ierr);
        k++;
        ierr = MatSetValues (A, 1, &j, 1, &k, &cc, INSERT_VALUES); CHKERRQ(ierr);
        k++;
        ierr = MatSetValues (A, 1, &j, 1, &k, &rc, INSERT_VALUES); CHKERRQ(ierr);
    }
    /* Left B.C. */
    j = 0;
    k = 0;
    ierr = MatSetValues (A, 1, &j, 1, &k, &cc, INSERT_VALUES); CHKERRQ(ierr);
    k++;
    ierr = MatSetValues (A, 1, &j, 1, &k, &rc, INSERT_VALUES); CHKERRQ(ierr);
    /* Right B.C. */
    j = n - 1;
    k = n - 2;
    ierr = MatSetValues (A, 1, &j, 1, &k, &lc, INSERT_VALUES); CHKERRQ(ierr);
    k++;
    ierr = MatSetValues (A, 1, &j, 1, &k, &cc, INSERT_VALUES); CHKERRQ(ierr);
    /* Finish matrix assembly */
    ierr = MatAssemblyEnd (A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatSetOption (A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE); CHKERRQ(ierr);
    /* Build r.h.s. */
    ierr = VecAssemblyBegin (Rhs); CHKERRQ(ierr);
    for (i = 1; i <= n; i++) {
        j = i - 1;
        t = T[i];
        ierr = VecSetValues (Rhs, 1, &j, &t, INSERT_VALUES); CHKERRQ(ierr);
    }
    ierr = VecCopy (Rhs, Sol); CHKERRQ(ierr);
    /* B.C. */
    ua = c*Ta;
    ub = c*Tb;
    /* Left */
    j = 0;
    ierr = VecSetValues (Rhs, 1, &j, &ua, ADD_VALUES); CHKERRQ(ierr);
    /* Right */
    j = n - 1;
    ierr = VecSetValues (Rhs, 1, &j, &ub, ADD_VALUES); CHKERRQ(ierr);
    ierr = VecAssemblyEnd (Rhs); CHKERRQ(ierr);

    /* Write the initial solution */
    sf_floatwrite (T, nx, solution);

    /* Loop in time */
    for (it = 1; it < nt; it++) {
        /* B.C. */
        T[0] = Ta;
        T[nx-1] = Tb;
        /* Solve the system */
        ierr = KSPSolve (Solver, Rhs, Sol); CHKERRQ(ierr);

        /* Move solution to r.h.s. */
        ierr = VecAssemblyBegin (Rhs); CHKERRQ(ierr);
        if (it != (nt - 1)) {
            ierr = VecCopy (Sol, Rhs); CHKERRQ(ierr);
            /* B.C. */
            j = 0;
            VecSetValues (Rhs, 1, &j, &ua, ADD_VALUES);
            j = n - 1;
            VecSetValues (Rhs, 1, &j, &ub, ADD_VALUES);
        }
        ierr = VecAssemblyEnd (Rhs); CHKERRQ(ierr);

        PetscFPrintf (comm, stderr, "Timestep #%d, t=%f\n", it, it*dt);

        /* Get solution from PETSc*/
        for (i = 1; i <= n; i++) {
            j = i - 1;
            VecGetValues (Sol, 1, &j, &t);
            T[i] = t;
        }
        /* Write the solution */
        sf_floatwrite (T, nx, solution);
    }

    /* Clean up */
    ierr = KSPDestroy (Solver); CHKERRQ(ierr);
    ierr = VecDestroy (Rhs); CHKERRQ(ierr);
    ierr = VecDestroy (Sol); CHKERRQ(ierr);
    ierr = PetscFinalize ();

    exit (0);
}
