/* Implicit solution of 2-D acoustic wave equation. */
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

#include "aimplfd2.h"

int main (int argc, char* argv[]) {
    int cpuid;
    int nz, nx, nt, it, ix, iz, niter;
    float oz, ox, ot, dz, dx, dt, sx, sz;
    float *u, *v, *f;
    bool fourth;
    sf_petsc_aimplfd2 aimplfd;
    /* I/O */
    sf_file usol=NULL, vel, src;
    /* PETSc */
    PetscErrorCode ierr;
	
    /* PETSc Initialization */
    ierr = PetscInitialize (&argc, &argv, 0, 0); CHKERRQ(ierr);
    MPI_Comm_rank (MPI_COMM_WORLD, &cpuid);
	
    sf_init (argc, argv);
	
    vel = sf_input ("in");
    /* Velocity */
    if (0 == cpuid)
        usol = sf_output ("out");
    /* Solution - wavefield */
	
    if (SF_FLOAT != sf_gettype (vel))
        sf_error ("Need float input");
    if (!sf_histint (vel, "n1", &nz)) sf_error ("No n1= in input");
    if (!sf_histfloat (vel, "d1", &dz)) sf_error ("No d1= in input");
    if (!sf_histfloat (vel, "o1", &oz)) oz = 0.;
    if (!sf_histint (vel, "n2", &nx)) sf_error ("No n2= in input");
    if (!sf_histfloat (vel, "d2", &dx)) sf_error ("No d2= in input");
    if (!sf_histfloat (vel, "o2", &ox)) ox = 0.;
	
    if (!sf_getstring ("src")) sf_error ("Need src=");
    /* Source wavelet */
    src = sf_input ("src");
    if (!sf_histint (src, "n1", &nt)) sf_error ("No n1= in src");
    if (!sf_histfloat (src, "d1", &dt)) sf_error ("No d1= in src");
    if (!sf_histfloat (src, "o1", &ot)) ot = 0.;
	
    if (!sf_histfloat (src, "o2", &sz)) sf_error ("No o2= in src");
    if (!sf_histfloat (src, "o3", &sx)) sf_error ("No o3= in src");
    iz = (sz - oz)/dz;
    ix = (sx - ox)/dx;
    if (iz < 0 || iz >= nz) sf_error ("Shot location is out of bounds");
    if (ix < 0 || ix >= nx) sf_error ("Shot location is out of bounds");
	
    if (!sf_getint ("niter", &niter)) niter = 10;
    /* Number of solver iterations */
	
    if (!sf_getbool ("fourth", &fourth)) fourth = true;
    /* Higher order flag */
	
    /* Velocity array */
    v = sf_floatalloc (nx*nz);
    /* Source array */
    f = sf_floatalloc (nt);
	
    /* Read velocity */
    sf_floatread (v, nx*nz, vel);
    /* Read source */
    sf_floatread (f, nt, src);
	
    /* Create time axis in output */
    if (0 == cpuid) {
        sf_putint (usol, "n3", nt);
        sf_putfloat (usol, "d3", dt);
        sf_putfloat (usol, "o3", ot);
    }
	
    PetscFPrintf (MPI_COMM_WORLD, stderr, "Initializing GMRES solver\n");
    
    aimplfd = sf_petsc_aimplfd2_init (nz, nx, dz, dx, dt, v, niter, fourth);
	
    free (v);
    /* Wavefield */
    u = sf_floatalloc (nx*nz);
	
    PetscFPrintf (MPI_COMM_WORLD, stderr, "Running GMRES solver\n");
	
    /* Loop in time */
    for (it = 0; it < nt; it++) {
        PetscFPrintf (MPI_COMM_WORLD, stderr, "Timestep #%d, t=%f\n", it, it*dt);
		
        sf_petsc_aimplfd2_next_step (aimplfd);
		
        sf_petsc_aimplfd2_add_source_ut1 (aimplfd, f[it], iz, ix);
        /* sf_petsc_aimplfd2_add_source_ut1 (aimplfd, f[it], iz, ix); */
        sf_petsc_aimplfd2_get_wavefield_ut2 (aimplfd, u);
		
        /* Write the solution */
        if (0 == cpuid)
            sf_floatwrite (u, nx*nz, usol);
    }
	
    /* Clean up */
    sf_petsc_aimplfd2_destroy (aimplfd);
    ierr = PetscFinalize ();
	
    exit (0);
}
