/* 2-D Tau tomography. */
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

#include <rsf.h>

static int nt, nx, nv, ntau;
float t0, dt;
static float **slopes; /* [ntau][nx] */
static float **tdata; /* [ntau][nx] */
static float **dtau; /* [ntau][nx] */
static float *v, *dv; /* [nt] */
static float *taus; /* [ntau] */

static void sf_tautomo_dtdtau () {
    int itau, ix, it;
    float t, tprev, p, q, td, tau, vel, dtdtau;

    /* Loop over events */
    for (itau = 0; itau < ntau; itau++) {
        /* Loop over offsets */
        for (ix = 0; ix < nx; ix++) {
            it = 0;
            t = 0;
            tprev = 0;
            tau = 0;
            p = slopes[itau][ix];
            td = tdata[itau][ix];
            dtdtau = 0;
            /* Loop until tau from data is reached */
            while (t < td && tprev < td) {
                tprev = t;
                vel = v[it];
                q = 1.0 - p*p*vel*vel;
                if (q < 0.0) break;
                t += dt/sqrt (q);
                /* Accept next step only if it is closer to tau from data */
                if ((t - td) < (td - tprev)) {
                    dtdtau = 1.0/sqrt (q);
                    tau += dt;
                    it++;
                }
            }
            dtau[itau][ix] = -dtdtau*(taus[itau] - tau);
        }
    }
}

static float sf_tautomo_dtdv (int itau, int ix, int iv) {
    if ((iv*dt - taus[itau]) > 0.5*dt)
        return 0.0;
    float p = slopes[itau][ix];
    float vel = v[iv];
    float q = 1 - p*p*vel*vel;
    if (q < 0.0)
        return 0.0;
    else
        return p*p*dt*vel/sqrt(q*q*q);
}

static void sf_tautomo_lop (bool adj, bool add, int nv, int ntau, float *dv, float *dtau)
{
    int iv, it, ix, itau;

    sf_adjnull (adj, add, nv, ntau, dv, dtau);

    for (iv = 0; iv < nv; iv++) {
        for (it = 0; it < ntau; it++) {
            ix = it % nx; /* Index of the offset */
            itau = it / nx; /* Index of the event */
            if (adj) {
                dv[iv] += sf_tautomo_dtdv (itau, ix, iv)*dtau[it];
            } else {
                dtau[it] += sf_tautomo_dtdv (itau, ix, iv)*dv[iv];
            }
        }
    }
}

static void sf_tautomo_rdata (float **data, float **pdata) {
    int tbase, itau, ix;
    float tfrac;
    for (itau = 0; itau < ntau; itau++) {
        tfrac = (taus[itau] - t0)/dt;
        tbase = (int)tfrac;
        tfrac -= tbase;
        for (ix = 0; ix < nx; ix++) {
            pdata[itau][ix] = data[ix][tbase]*(1 - tfrac) + data[ix][tbase + 1]*tfrac;
        }
    }
}

int main (int argc, char* argv[]) {
    bool verb = false;
    int i, k, l, iter, npicks, niter = 3;
    float **data;
    sf_file in, out, time, vel0, picks;

    sf_init (argc, argv);
    in = sf_input ("in");
    /* Flattened slopes */
    out = sf_output ("out");
    /* Estimated velocities */

    if (!sf_histint (in, "n1", &nt)) sf_error ("Need n1= in in");
    if (!sf_histint (in, "n2", &nx)) sf_error ("Need n2= in in");
    if (SF_FLOAT != sf_gettype (in)) sf_error ("Need float input");
    if (!sf_histfloat (in, "d1", &dt)) sf_error ("Need d1= in in");
    if (!sf_histfloat (in, "o1", &t0)) t0 = 0.0;
    nv = nt;

    if (NULL == sf_getstring ("time")) sf_error ("Need time= input");
    /* Flattened traveltimes */
    time = sf_input ("time");
    if (NULL == sf_getstring ("vel0")) sf_error ("Need vel0= input");
    /* Initial velocities */
    vel0 = sf_input ("vel0");
    if (NULL == sf_getstring ("picks")) sf_error ("Need picks= input");
    /* Key event picks */
    picks = sf_input ("picks");
    if (!sf_histint (picks, "n1", &npicks)) sf_error ("Need n1= in picks=");
    
    if (!sf_getint ("niter", &niter)) niter = 3;
    /* Number of iterations */
    if (!sf_getbool ("verb", &verb)) verb = false;
    /* Verbosity flag */

    sf_unshiftdim (in, out, 2);
    /* Collapse offset dimension */

    data = sf_floatalloc2 (nt, nx);
    taus = sf_floatalloc (nt);
    v = sf_floatalloc (nt);
    dv = sf_floatalloc (nt);
    for (i = 0; i < nt; i++)
        taus[i] = SF_HUGE;

    k = 0; i = 0;
    sf_floatread (&taus[0], 1, picks);
    while (i < npicks) { /* Loop over picks */
        k++;
        /* If end of picks for current ensemble or whole dataset */
        if ((npicks - 1) == i ||  taus[k] < taus[k - 1]) {
            if (verb) sf_warning ("Processing CMP #%d", i + 1);
            ntau = k;
            slopes = sf_floatalloc2 (nx, ntau);
            tdata = sf_floatalloc2 (nx, ntau);
            dtau = sf_floatalloc2 (nx, ntau);
            /* Read and extract slopes according to picks */
            sf_floatread (data[0], nt*nx, in);
            sf_tautomo_rdata (data, slopes);
            /* Read and extract times according to picks */
            sf_floatread (data[0], nt*nx, time);
            sf_tautomo_rdata (data, tdata);
            /* Read initial velocities */
            sf_floatread (v, nt, vel0);
            /* Dot test */
/*
            double dot1[2], dot2[2];
            if (verb) sf_warning ("Dot test");
            sf_dot_test (sf_tautomo_lop, nv, ntau*nx, dot1, dot2);
            sf_warning ("%12.8f ? %12.8f\n", dot1[0], dot1[1]);
            sf_warning ("%12.8f ? %12.8f\n", dot2[0], dot2[1]);
*/
            /* Prepare right-hand side */
            for (iter = 0; iter < niter; iter++) {
                memset (dv, 0, nt*sizeof(float));
                if (verb) sf_warning ("Linear iteration %d", iter + 1);
                sf_tautomo_dtdtau ();
                /* Solve the linear system */
                sf_solver (sf_tautomo_lop, sf_cgstep, nv, ntau*nx, dv, dtau[0], niter,
                           "verb", verb, "end");
                for (l = 0; l < nt; l++)
                    v[l] += dv[l];
            }
            /* Write results */
            sf_floatwrite (v, nt, out);
            /* Move on to the next ensemble */
            taus[0] = taus[k];
            taus[1] = SF_HUGE;
            k = 0;
            free (slopes[0]);
            free (tdata[0]);
            free (dtau[0]);
            free (slopes);
            free (tdata);
            free (dtau);
        } else /* Read next pick */
            sf_floatread (&taus[k], 1, picks);
        i++;
    }

    sf_fileclose (picks);
    sf_fileclose (vel0);
    sf_fileclose (time);

    return 0;
}
