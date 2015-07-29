/* Implicit solution of 1-D acoustic wave equation. */
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

#include "aimplfd1.h"

int main (int argc, char* argv[]) {
    int cpuid;
    int nx, nt, it, is, niter;
    float ox, ot, dx, dt, sx;
    float *u, *v, *f;
    sf_petsc_aimplfd1 aimplfd;
    /* I/O */
    sf_file usol=NULL, vel, src;
    /* PETSc */
    PetscErrorCode ierr;

    /* PETSc Initialization */
    ierr = PetscInitialize (&argc, &argv, 0, 0); CHKERRQ(ierr);
    MPI_Comm_rank (MPI_COMM_WORLD, &cpuid);
    sf_warning ("CPU id %d", cpuid);

    sf_init (argc, argv);

    vel = sf_input ("in");
    /* Velocity */
    if (0 == cpuid)
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

    if (!sf_getint ("niter", &niter)) niter = 10;
    /* Number of solver iterations */

    /* Velocity array */
    v = sf_floatalloc (nx);
    /* Source array */
    f = sf_floatalloc (nt);

    /* Read velocity */
    sf_floatread (v, nx, vel);
    /* Read source */
    sf_floatread (f, nt, src);

    /* Create time axis in output */
    if (0 == cpuid) {
        sf_putint (usol, "n2", nt);
        sf_putfloat (usol, "d2", dt);
        sf_putfloat (usol, "o2", ot);
    }

    PetscFPrintf (MPI_COMM_WORLD, stderr, "Initializing GMRES solver\n");
    
    aimplfd = sf_petsc_aimplfd1_init (nx, dx, dt, v, niter);

    free (v);
    /* Wavefield */
    u = sf_floatalloc (nx);

    PetscFPrintf (MPI_COMM_WORLD, stderr, "Running GMRES solver\n");

    /* Loop in time */
    for (it = 0; it < nt; it++) {
        PetscFPrintf (MPI_COMM_WORLD, stderr, "Timestep #%d, t=%f\n", it, it*dt);

        sf_petsc_aimplfd1_next_step (aimplfd);
        sf_petsc_aimplfd1_add_source (aimplfd, f[it], is);
        sf_petsc_aimplfd1_get_wavefield (aimplfd, u);

        /* Write the solution */
        if (0 == cpuid)
            sf_floatwrite (u, nx, usol);
    }

    /* Clean up */
    sf_petsc_aimplfd1_destroy (aimplfd);
    ierr = PetscFinalize ();

    exit (0);
}
