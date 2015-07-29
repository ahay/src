/* Implicit solution of 2-D acoustic wave equation, compatibility interface with sfawefd2d */
/*
  Copyright (C) 2007 Colorado School of Mines
  
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
#include "fdutil.h"

int main (int argc, char* argv[])
{
    bool verb, fsrf, snap, expl, dabc; 
    int  jsnap, ntsnap, jdata;

    /* I/O files */
    sf_file Fwav = NULL; /* wavelet   */
    sf_file Fsou = NULL; /* sources   */
    sf_file Frec = NULL; /* receivers */
    sf_file Fvel = NULL; /* velocity  */
    sf_file Fden = NULL; /* density   */
    sf_file Fdat = NULL; /* data      */
    sf_file Fwfl = NULL; /* wavefield */

    /* cube axes */
    sf_axis at, az, ax;
    sf_axis as, ar;

    int cpuid;
    int nt, nz, nx, ns, nr;
    int it, ia;
    float dt, dz, dx;

    /* FDM structure */
    /* Wee keep these structures for compatibility,
       as soon as PETSC version is finished,
       this will be removed */
    fdm2d fdm = NULL;

    /* I/O arrays */
    float *ww = NULL;           /* wavelet   */
    pt2d  *ss = NULL;           /* sources   */
    pt2d  *rr = NULL;           /* receivers */
    float *dd = NULL;           /* data      */

    float **u,**v;

    /* linear interpolation weights/indices */
    lint2d cs,cr;

    /* wavefield cut params */
    sf_axis acz = NULL, acx = NULL;
    int   nqz, nqx;
    float oqz, oqx;
    float dqz, dqx;
    float **uc = NULL;

    sf_petsc_aimplfd2 aimplfd;
    PetscErrorCode ierr;
    /* PETSc Initialization */
    ierr = PetscInitialize (&argc, &argv, 0, 0); CHKERRQ(ierr);
    MPI_Comm_rank (MPI_COMM_WORLD, &cpuid);

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init (argc, argv);

    /*------------------------------------------------------------*/

    if (!sf_getbool ("verb", &verb)) verb = false; /* verbosity flag */
    if (!sf_getbool ("snap", &snap)) snap = false; /* wavefield snapshots flag */
    if (!sf_getbool ("free", &fsrf)) fsrf = false; /* free surface flag */
    if (!sf_getbool ("expl", &expl)) expl = false; /* "exploding reflector" */
    if (!sf_getbool ("dabc", &dabc)) dabc = false; /* absorbing BC */
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* I/O files */
    Fwav = sf_input ("in" ); /* wavelet   */
    Fvel = sf_input ("vel"); /* velocity  */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */
    Fwfl = sf_output("wfl"); /* wavefield */
    Fdat = sf_output("out"); /* data      */
    Fden = sf_input ("den"); /* density   */

    /*------------------------------------------------------------*/
    /* axes */
    at = sf_iaxa (Fwav,2); sf_setlabel (at,"t"); if (verb && 0 == cpuid) sf_raxa (at); /* time */
    az = sf_iaxa (Fvel,1); sf_setlabel (az,"z"); if (verb && 0 == cpuid) sf_raxa (az); /* depth */
    ax = sf_iaxa (Fvel,2); sf_setlabel (ax,"x"); if (verb && 0 == cpuid) sf_raxa (ax); /* space */

    as = sf_iaxa (Fsou, 2); sf_setlabel (as, "s"); if (verb && 0 == cpuid) sf_raxa (as); /* sources */
    ar = sf_iaxa (Frec, 2); sf_setlabel (ar, "r"); if (verb && 0 == cpuid) sf_raxa (ar); /* receivers */

    nt = sf_n (at); dt = sf_d (at);
    nz = sf_n (az); dz = sf_d (az);
    nx = sf_n (ax); dx = sf_d (ax);

    ns = sf_n (as);
    nr = sf_n (ar);
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* other execution parameters */
    if (!sf_getint ("jdata", &jdata)) jdata = 1;
    if (snap) {  /* save wavefield every *jsnap* time steps */
        if (!sf_getint ("jsnap", &jsnap)) jsnap = nt;
    }
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    /*if (!sf_getint("nb",&nb) || nb<NOP) nb=NOP;*/
    fdm = fdutil_init (verb, true, az, ax, 0, 1);

    /*------------------------------------------------------------*/

    if (0 == cpuid) {
        sf_setn (az, fdm->nzpad); sf_seto (az, fdm->ozpad); if(verb) sf_raxa (az);
        sf_setn (ax, fdm->nxpad); sf_seto (ax, fdm->oxpad); if(verb) sf_raxa (ax);
        /*------------------------------------------------------------*/

        /*------------------------------------------------------------*/
        /* setup output data header */
        sf_oaxa (Fdat, ar, 1);

        sf_setn (at, nt/jdata);
        sf_setd (at, dt*jdata);
        sf_oaxa (Fdat, at, 2);
    }

    /* setup output wavefield header */
    if (snap && 0 == cpuid) {
        if (!sf_getint ("nqz", &nqz)) nqz = sf_n (az);
        if (!sf_getint ("nqx", &nqx)) nqx = sf_n (ax);

        if (!sf_getfloat ("oqz", &oqz)) oqz = sf_o (az);
        if (!sf_getfloat ("oqx", &oqx)) oqx = sf_o (ax);

        dqz = sf_d (az);
        dqx = sf_d (ax);

        acz = sf_maxa (nqz, oqz, dqz); sf_raxa (acz);
        acx = sf_maxa (nqx, oqx, dqx); sf_raxa (acx);
        /* check if the imaging window fits in the wavefield domain */

        uc = sf_floatalloc2 (sf_n (acz), sf_n (acx));

        ntsnap = 0;
        for (it = 0; it < nt; it++) {
            if (it % jsnap == 0)
                ntsnap++;
        }
        sf_setn (at, ntsnap);
        sf_setd (at, dt*jsnap);
        if (verb)
            sf_raxa(at);

        sf_oaxa (Fwfl, acz, 1);
        sf_oaxa (Fwfl, acx, 2);
        sf_oaxa (Fwfl, at,  3);
    }

    if (expl) {
        ww = sf_floatalloc (1);
    } else {
        ww = sf_floatalloc (ns);
    }
    dd = sf_floatalloc (nr);

    /*------------------------------------------------------------*/
    /* setup source/receiver coordinates */
    ss = (pt2d*) sf_alloc (ns, sizeof (*ss));
    rr = (pt2d*) sf_alloc (nr, sizeof (*rr));

    pt2dread1 (Fsou, ss, ns, 2); /* read (x,z) coordinates */
    pt2dread1 (Frec, rr, nr, 2); /* read (x,z) coordinates */

    cs = lint2d_make (ns, ss, fdm);
    cr = lint2d_make (nr, rr, fdm);

    v = sf_floatalloc2 (nz, nx);
    u = sf_floatalloc2 (nz, nx);

    /* input velocity */
    sf_floatread (v[0], nz*nx, Fvel);

    /*------------------------------------------------------------*/
    PetscFPrintf (MPI_COMM_WORLD, stderr, "Initializing GMRES solver\n");   
    aimplfd = sf_petsc_aimplfd2_init (nz, nx, dz, dx, dt, &v[0][0], 100, true);

    /*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
    for (it = 0; it < nt; it++) {
        PetscFPrintf (MPI_COMM_WORLD, stderr, "Timestep #%d, t=%f\n", it, it*dt);

        sf_petsc_aimplfd2_next_step (aimplfd);

        /* inject acceleration source */
        if (expl) {
            sf_floatread (ww, 1, Fwav);
            for (ia = 0; ia < cs->n; ia++) {
                sf_petsc_aimplfd2_add_source_ut1 (aimplfd, ww[0], cs->jz[ia], cs->jx[ia]);
            }
        } else {
            sf_floatread (ww, ns, Fwav);
            for (ia = 0; ia < cs->n; ia++) {
            /*
                PetscFPrintf (MPI_COMM_WORLD, stderr, "Source #%d [%d, %d], f=%f\n", ia, cs->jz[ia], cs->jx[ia], ww[0]);
                */
                sf_petsc_aimplfd2_add_source_ut1 (aimplfd, ww[ia], cs->jz[ia], cs->jx[ia]);
            }
        }
        sf_petsc_aimplfd2_get_wavefield_ut2 (aimplfd, &u[0][0]);

        /* extract data */
        /*
        lint2d_extract (u, dd, cr);
        */
        for (ia = 0; ia < cr->n; ia++) {
            dd[ia] = u[cr->jx[ia]][cr->jz[ia]];
        }

        if (snap && it % jsnap == 0 && 0 == cpuid) {
            cut2d (u, uc, fdm, acz, acx);
            sf_floatwrite (uc[0], sf_n(acz)*sf_n(acx), Fwfl);
        }
        if (it % jdata == 0 && 0 == cpuid)
            sf_floatwrite (dd, nr, Fdat);
    }

    exit (0);
}

