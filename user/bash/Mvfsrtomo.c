/* Global 2-D tomography with VFSR. */
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

#include <math.h>
#include <rsf.h>

#include "eikfswp.h"

/* Interpolates 2D velocity to a bigger grid */
static void sf_vel_spline_interp_2d (double *v0, float *v,
                                     int n01, int n02,
                                     int n1, int n2) {
    int i, j;
    double val;

    Ugrid y_grid, z_grid;
    y_grid.start = 0.0; y_grid.end = n2 - 1 + 0.0001; y_grid.num = n02;
    z_grid.start = 0.0; z_grid.end = n1 - 1 + 0.0001; z_grid.num = n01;
    BCtype_d yBC, zBC;
    yBC.lCode = yBC.rCode = NATURAL;
    zBC.lCode = zBC.rCode = NATURAL;

    UBspline_2d_d *spline = create_UBspline_2d_d (y_grid, z_grid, yBC, zBC, v0);

    for (j = 0; j < n2; j++) {
        for (i = 0; i < n1; i++) {
            eval_UBspline_2d_d (spline, (double)j, (double)i, &val);
            v[j * n1 + i] = val;
        }
    }

    destroy_Bspline (spline);
}

typedef enum { TOP, BOTTOM, LEFT, RIGHT } tt_pos;

typedef struct {
    float **s; /* Shot coordinates */
    int     ns; /* Number of shots */
    float  *v; /* Storage for interpolated velocity */
    float  *t; /* Storage for calculated traveltime */
    float **tobs; /* Observed traveltimes */
    tt_pos *ttime_pos; /* Traveltime location (top/bottom/etc) */
    float   d1, d2;
    float   o1, o2;
    int     n1, n2;
    int     np1, np2; /* Number of reduced samples */
    int     niter; /* Number of sweeping iterations */
    bool    done; /* if TRUE, do not evaluate cost function */
    float   cost_eps; /* Minimum value for the cost function */
} sf_vfsr_tomo_data;

/* Cost function to be optimized -
   it is an L2-norm between the observed traveltimes
   and calculated traveltimes for all shots */
double sf_vfsr_tomo_cost_func (double *cost_parameters,
                               double *parameter_lower_bound,
                               double *parameter_upper_bound,
                               int *cost_flag,
                               VFSR_DEFINES *USER_OPTIONS,
                               void *user_data) {
    sf_vfsr_tomo_data *data = (sf_vfsr_tomo_data*)user_data;
    int i, j, num = 1, start = 0, inc = 1;
    double diff = 0.0;

    if (data->done) {
        *cost_flag = FALSE;
        return 0;
    }

    sf_vel_spline_interp_2d (cost_parameters, data->v,
                             data->np1, data->np2,
                             data->n1, data->n2);

    /* loop over shots */
    for (i = 0; i < data->ns; i++) {
        if (false == sf_init_fast_sweep (data->t,
                                         1, data->n2, data->n1,
                                         0, data->o2, data->o1,
                                         1, data->d2, data->d1,
                                         data->s[i][0], data->s[i][1], 0))
            sf_error ("Incorrect shot location");

        sf_run_fast_sweep (data->t, data->v, data->niter,
                           1, data->n2, data->n1,
                           0, data->o2, data->o1,
                           0, data->d2, data->d1);
        switch (data->ttime_pos[i]) {
            case TOP:
                num = data->n2;
                start = 0;
                inc = data->n1;
                break;
            case BOTTOM:
                num = data->n2;
                start = data->n1 - 1;
                inc = data->n1;
                break;
            case LEFT:
                num = data->n1;
                start = 0;
                inc = 1;
                break;
            case RIGHT:
                num = data->n1;
                start = (data->n2 - 1) * data->n1;
                inc = 1;
                break;
        }
        for (j = 0; j < num; j++) {
            diff += pow (data->t[start + j * inc] -
                         data->tobs[i][j], 2.0);
        }
    }

    diff = sqrt (diff);

    if (diff < data->cost_eps)
        data->done = true;

    *cost_flag = TRUE;
    return diff;
}

int main (int argc, char* argv[]) {
    int n1, n2, nt, i, j;
    int nshot, niter, np1, np2, dp1, dp2, np;
    float o1, o2, d1, d2, vmin, vmax, eps;
    float **s, **tobs = NULL, *t, *v, *vtmp;
    tt_pos *ttime_pos = NULL;
    bool isvel;

    /* I/O */
    char *sfile, *tfile, *lfile;
    sf_file vel0, vel, time, shots, locs;

    /* VFSR data */
    int exit_status, *param_type;
    double *params, *param_min, *param_max, *param_tans, *param_curv;
    VFSR_DEFINES *vfsr_user_options;
    sf_vfsr_tomo_data vfsr_user_data;

    sf_init (argc, argv);
    vel0 = sf_input ("in");
    vel = sf_output ("out");

    if (SF_FLOAT != sf_gettype (vel0))
        sf_error("Need float input");
    if (!sf_histint (vel0, "n1", &n1)) sf_error ("No n1= in input");
    if (!sf_histint (vel0, "n2", &n2)) sf_error ("No n2= in input");

    if (!sf_histfloat (vel0, "d1", &d1)) sf_error ("No d1= in input");
    if (!sf_histfloat (vel0, "d2", &d2)) sf_error ("No d2= in input");

    if (!sf_histfloat (vel0, "o1", &o1)) o1 = 0.;
    if (!sf_histfloat (vel0, "o2", &o2)) o2 = 0.;

    if (!sf_getbool ("vel", &isvel)) isvel = true;
    /* if y, the input is velocity; n - slowness */

    if (!sf_getint ("niter", &niter)) niter = 2;
    /* number of sweeping iterations */

    if (!sf_getfloat ("eps", &eps)) sf_error ("Need eps=");
    /* minimum value for the cost function to be reached */

    if (!sf_getfloat ("vmin", &vmin)) sf_error ("Need vmin=");
    /* minimum velocity */
    if (!sf_getfloat ("vmax", &vmax)) sf_error ("Need vmax=");
    /* maximum velocity */

    if (!sf_getint ("np1", &np1)) sf_error ("Need np1=");
    /* number of skipped velocity values in n1 direction */
    if (((n1 - 1) / np1) * np1 != (n1 - 1))
        sf_error ("np1 does not fit n1");
    dp1 = np1;
    np1 = (n1 - 1) / np1 + 1;
    if (!sf_getint ("np2", &np2)) sf_error ("Need np2=");
    /* number of skipped velocity values in n2 direction */
    if (((n1 - 1) / np2) * np2 != (n2 - 1))
        sf_error ("np2 does not fit n2");
    dp2 = np2;
    np2 = (n2 - 1) / np2 + 1;
    np = np1 * np2;

    sfile = sf_getstring ("shotfile");
    /* File with shot locations (n2=number of shots, n1=2) */

    if (NULL != sfile) {
        shots = sf_input ("shotfile");

        if (SF_FLOAT != sf_gettype (shots)) 
            sf_error ("Need float shotfile");
        if (!sf_histint (shots, "n2", &nshot))
            sf_error ("No n2= in shotfile");
        if (!sf_histint (shots, "n1", &i) || i < 2)
            sf_error ("Need n1>=2 in shotfile");

        s = sf_floatalloc2 (i, nshot);
        sf_floatread (s[0], nshot * i, shots);

        free (sfile);
    } else {
        nshot = 1;

        s = sf_floatalloc2 (2, nshot);

        if (!sf_getfloat ("zshot", &s[0][0])) s[0][0] = 0.; 
        /* Shot location (used if no shotfile) */
        if (!sf_getfloat ("yshot", &s[0][1])) s[0][1] = o2 + 0.5*(n2 - 1)*d2;
    }

    tfile = sf_getstring ("timefile");
    /* File with observed traveltimes (n2=number of shots) */

    if (NULL != tfile) {
        time = sf_input ("timefile");

        if (SF_FLOAT != sf_gettype (time)) 
            sf_error ("Need float timefile");
        if (!sf_histint (time, "n1", &nt))
            sf_error ("No n1= in timefile");
        if (nt != n1 && nt != n2)
            sf_error ("Need n1 in timefile to be equal to n1 or n2 of input");
        if (!sf_histint (time, "n2", &i) || i != nshot)
            sf_error ("Need n2 in timefile to be equal to the number of shots");

        tobs = sf_floatalloc2 (nt, nshot);
        sf_floatread (tobs[0], nt * nshot, time);
    } else
        sf_error ("Need timefile=");

    lfile = sf_getstring ("locfile");
    /* File with encoded locations of the traveltimes (n2=number of shots, n1=1),
       encoding: 0 - top, 1 - bottom, 2 - left, 3 - right */

    if (NULL != lfile) {
        locs = sf_input ("locfile");

        if (SF_INT != sf_gettype (locs)) 
            sf_error ("Need integer locfile");
        if (!sf_histint (locs, "n1", &i))
            sf_error ("No n1= in locfile");
        if (i != 1)
            sf_error ("Need n1=1 in locfile");
        if (!sf_histint (locs, "n2", &i) || i != nshot)
            sf_error ("Need n2 in locfile to be equal to number of shots");

        ttime_pos = (tt_pos*)sf_intalloc (nshot);
        sf_intread ((int*)ttime_pos, nshot, locs);
    } else
        sf_error ("Need locfile=");

    t  = sf_floatalloc (n1 * n2);
    v  = sf_floatalloc (n1 * n2);
    vtmp = sf_floatalloc (n1 * n2);

    sf_floatread (v, n1 * n2, vel0);
    /* transform velocity to slowness */
    if (isvel) {
        for (i = 0; i < n1 * n2; i++) {
            v[i] = 1. / v[i];
        }
        t[0] = 1.0 / vmax;
        t[1] = 1.0 / vmin;
        vmin = t[0];
        vmax = t[1];
    }

    param_type = sf_intalloc (np);
    params = sf_alloc (np, sizeof (double));
    param_min = sf_alloc (np, sizeof (double));
    param_max = sf_alloc (np, sizeof (double));
    param_tans = sf_alloc (np, sizeof (double));
    param_curv = sf_alloc (np * np, sizeof (double));

    sf_warning ("Min, max: %g %g", vmin, vmax);
    sf_warning ("Using %dx%d samples", np1, np2);

    for (i = 0; i < np; i++) {
        param_min[i] = vmin;
        param_max[i] = vmax;
        param_type[i] = REAL_TYPE;
    }

    for (j = 0; j < np2; j++) {
        for (i = 0; i < np1; i++) {
            params[j * np1 + i] = v[j * dp2 * n1 + i * dp1];
        }
    }

    vfsr_user_data.s = s;
    vfsr_user_data.ns = nshot;
    vfsr_user_data.v = vtmp;
    vfsr_user_data.t = t;
    vfsr_user_data.tobs = tobs;
    vfsr_user_data.ttime_pos = ttime_pos;
    vfsr_user_data.d1 = d1;
    vfsr_user_data.d2 = d2;
    vfsr_user_data.o1 = o1;
    vfsr_user_data.o2 = o2;
    vfsr_user_data.n1 = n1;
    vfsr_user_data.n2 = n2;
    vfsr_user_data.np1 = np1;
    vfsr_user_data.np2 = np2;
    vfsr_user_data.niter = niter;
    vfsr_user_data.done = false;
    vfsr_user_data.cost_eps = eps;

    vfsr_user_options = vfsr_get_new_defines ();
    vfsr_user_options->USER_INITIAL_PARAMETERS = TRUE;

    vfsr_std_rng (sf_vfsr_tomo_cost_func,
                  np, param_type, params,
                  param_min, param_max,
                  param_tans, param_curv,
                  &exit_status,
                  vfsr_user_options,
                  (void*)&vfsr_user_data);

    sf_warning ("VFSR exit code: %d", exit_status);

    sf_vel_spline_interp_2d (params, v,
                             np1, np2,
                             n1, n2);

    if (isvel) {
        for (i = 0; i < n1 * n2; i++) {
            v[i] = 1. / v[i];
        }
    }

    sf_floatwrite (v, n1 * n2, vel);

    exit (0);
}
