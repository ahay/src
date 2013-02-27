/* 5-D phase-space escape values grid consisting of supercells */
/*
  Copyright (C) 2012 University of Texas at Austin

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

#include <errno.h>
#include <sys/mman.h>

#include <rsf.h>

#ifndef _esc_scgrid3_h

#include "esc_tracer3.h"
/*^*/

typedef struct EscSCgrid3 *sf_esc_scgrid3;
/* abstract data type */
/*^*/

#endif

#include "einspline.h"

struct EscSCgrid3 {
    size_t               n, offs, is, ir;
    int                  nz, nx, ny, na, nb;
    float                oz, ox, oy, oa, ob;
    float                dz, dx, dy, da, db;
    float                zmin, zmax, xmin, xmax, ymin, ymax;
    bool                 quad;
    multi_UBspline_3d_s *scsplines;
    unsigned char       *mmaped;
    sf_esc_point3        esc_point;
    sf_esc_tracer3       esc_tracer;
};
/* concrete data type */

sf_esc_scgrid3 sf_esc_scgrid3_init (sf_file scgrid, sf_esc_tracer3 esc_tracer, bool verb)
/*< Initialize object >*/
{
    size_t nc, nnc = 0;
    int i, ia, ib;
    size_t offs;
    FILE *stream;
    sf_esc_scgrid3 esc_scgrid = (sf_esc_scgrid3)sf_alloc (1, sizeof (struct EscSCgrid3));

    if (sf_gettype (scgrid) != SF_UCHAR)
        sf_error ("Wrong data type in supercell file");

    if (!sf_histlargeint (scgrid, "n1", (off_t*)&esc_scgrid->n)) sf_error ("No n1= in supercell file");

    if (!sf_histint (scgrid, "Nz", &esc_scgrid->nz)) sf_error ("No Nz= in supercell file");
    if (!sf_histfloat (scgrid, "Dz", &esc_scgrid->dz)) sf_error ("No Dz= in supercell file");
    if (!sf_histfloat (scgrid, "Oz", &esc_scgrid->oz)) sf_error ("No Oz= in supercell file");
    if (!sf_histint (scgrid, "Nx", &esc_scgrid->nx)) sf_error ("No Nx= in supercell file");
    if (!sf_histfloat (scgrid, "Dx", &esc_scgrid->dx)) sf_error ("No Dx= in supercell file");
    if (!sf_histfloat (scgrid, "Ox", &esc_scgrid->ox)) sf_error ("No Ox= in supercell file");
    if (!sf_histint (scgrid, "Ny", &esc_scgrid->ny)) sf_error ("No Ny= in supercell file");
    if (!sf_histfloat (scgrid, "Dy", &esc_scgrid->dy)) sf_error ("No Dy= in supercell file");
    if (!sf_histfloat (scgrid, "Oy", &esc_scgrid->oy)) sf_error ("No Oy= in supercell file");
    if (!sf_histint (scgrid, "n2", &esc_scgrid->nb)) sf_error ("No n2= in supercell file");
    if (!sf_histfloat (scgrid, "d2", &esc_scgrid->db)) sf_error ("No d2= in supercell file");
    if (!sf_histfloat (scgrid, "o2", &esc_scgrid->ob)) sf_error ("No o2= in supercell file");
    if (!sf_histint (scgrid, "n3", &esc_scgrid->na)) sf_error ("No n3= in supercell file");
    if (!sf_histfloat (scgrid, "d3", &esc_scgrid->da)) sf_error ("No d3= in supercell file");
    if (!sf_histfloat (scgrid, "o3", &esc_scgrid->oa)) sf_error ("No o3= in supercell file");

    esc_scgrid->zmin = esc_scgrid->oz;
    esc_scgrid->zmax = esc_scgrid->oz + (esc_scgrid->nz - 1)*esc_scgrid->dz;
    esc_scgrid->xmin = esc_scgrid->ox;
    esc_scgrid->xmax = esc_scgrid->ox + (esc_scgrid->nx - 1)*esc_scgrid->dx;
    esc_scgrid->ymin = esc_scgrid->oy;
    esc_scgrid->ymax = esc_scgrid->oy + (esc_scgrid->ny - 1)*esc_scgrid->dy;

    stream = sf_filestream (scgrid);
    if (stdin == stream)
        sf_error ("Can not mmap supercell file from stdin");

    esc_scgrid->offs = ftello (stream);

    esc_scgrid->scsplines = (multi_UBspline_3d_s*)sf_alloc ((size_t)esc_scgrid->na*(size_t)esc_scgrid->nb,
                                                            sizeof(multi_UBspline_3d_s));
    esc_scgrid->mmaped = (unsigned char*)mmap (NULL, (size_t)esc_scgrid->offs +
                                                     (size_t)esc_scgrid->n*
                                                     (size_t)esc_scgrid->na*
                                                     (size_t)esc_scgrid->nb,
                                               PROT_READ, MAP_SHARED,
                                               fileno (stream), 0);
    nc = esc_scgrid->offs;
    for (ia = 0; ia < esc_scgrid->na; ia++) {
        for (ib = 0; ib < esc_scgrid->nb; ib++) {
            /* Copy spline structure to a separate place to
               make it modifiable */
            memcpy ((void*)&esc_scgrid->scsplines[ia*esc_scgrid->nb + ib],
                    (void*)&esc_scgrid->mmaped[nc], sizeof(multi_UBspline_3d_s));
            nc += sizeof(multi_UBspline_3d_s);
            /* Initialize pointer to coefficients with a correct
               new address in the mmaped area */
            esc_scgrid->scsplines[ia*esc_scgrid->nb + ib].coefs = (float*)&esc_scgrid->mmaped[nc];
            nc += esc_scgrid->scsplines[ia*esc_scgrid->nb + ib].nc;
            nnc += esc_scgrid->scsplines[ia*esc_scgrid->nb + ib].nc;
        }
    }

    if (verb) {
        sf_warning ("Spatial domain dimensions: nz=%d, z=[%g, %g]", esc_scgrid->nz,
                    esc_scgrid->zmin, esc_scgrid->zmax);
        sf_warning ("Spatial domain dimensions: nx=%d, x=[%g, %g]", esc_scgrid->nx,
                    esc_scgrid->xmin, esc_scgrid->xmax);
        sf_warning ("Spatial domain dimensions: ny=%d, y=[%g, %g]", esc_scgrid->ny,
                    esc_scgrid->ymin, esc_scgrid->ymax);
        sf_warning ("Angular domain dimensions: nb=%d, b=[%g, %g]", esc_scgrid->nb,
                    esc_scgrid->ob, esc_scgrid->ob + (esc_scgrid->nb - 1)*esc_scgrid->db);
        sf_warning ("Angular domain dimensions: na=%d, a=[%g, %g]", esc_scgrid->na,
                    esc_scgrid->oa, esc_scgrid->oa + (esc_scgrid->na - 1)*esc_scgrid->da);
        sf_warning ("%lu supercells read, %g Mb of spline coefficients accumulated",
                    (size_t)esc_scgrid->na*(size_t)esc_scgrid->nb, nnc*1e-6);
    }

    esc_scgrid->da *= SF_PI/180.0;
    esc_scgrid->oa *= SF_PI/180.0;
    esc_scgrid->db *= SF_PI/180.0;
    esc_scgrid->ob *= SF_PI/180.0;

    esc_scgrid->is = 0; /* Total number of interpolation steps */
    esc_scgrid->ir = 0; /* Total number of processed points */

    esc_scgrid->esc_tracer = esc_tracer;
    esc_scgrid->esc_point = sf_esc_point3_init ();

    return esc_scgrid;
}

void sf_esc_scgrid3_close (sf_esc_scgrid3 esc_scgrid, bool verb)
/*< Destroy object >*/
{
    if (verb)
        sf_warning ("%lu points processed, %g interpolation steps per point performed",
                    esc_scgrid->ir, (float)esc_scgrid->is/(float)esc_scgrid->ir);
    munmap (esc_scgrid->mmaped, (size_t)esc_scgrid->offs +
                                (size_t)esc_scgrid->n*
                                (size_t)esc_scgrid->na*
                                (size_t)esc_scgrid->nb);
    free (esc_scgrid->scsplines);
    sf_esc_point3_close (esc_scgrid->esc_point);
    free (esc_scgrid);
}

static int sf_cram_scgrid3_ftoi (float v, float f0, float df, float *f) {
    int i;
    float ff;

    ff = (v - f0)/df;
    i = floorf (ff);
    if (f)
        *f = ff - (float)i;
    return i;
}

/* Convert normalized phase vector components into polar angles */
static void sf_esc_scgrid3_p_to_ab (sf_esc_scgrid3 esc_scgrid,
                                    float pz, float px, float py,
                                    float *a, float *b) {
    if (pz < -1.0) pz = -1.0;
    if (pz > 1.0) pz = 1.0;
    if (px < -1.0) px = -1.0;
    if (px > 1.0) px = 1.0;
    if (py < -1.0) py = -1.0;
    if (py > 1.0) py = 1.0;

    *b = acosf (pz);
    *a = acosf (px);
    if (py < 0.0)
        *a = 2.0*SF_PI - *a;

    /* Make sure that a is in [0;2PI] and b is in [0;PI] */
    if (*a < 0.0)
        *a += 2.0*SF_PI;
    else if (*a > 2.0*SF_PI)
        *a -= 2.0*SF_PI;
    if (*b < 0.0) {
        *a += SF_PI;
        *b = fabsf (*b);
    } else if (*b > SF_PI) {
        *a -= SF_PI;
        *b = 2.0*SF_PI - *b;
    }  
    if (*a < 0.0)
        *a += 2.0*SF_PI;
    else if (*a > 2.0*SF_PI)
        *a -= 2.0*SF_PI;
}

/* Flip azimuth by +/-PI (this happens when inclination is moved back to
   its interval of [0;PI) */
static int sf_esc_scgrid3_flip_ia (sf_esc_scgrid3 esc_scgrid, int ia) {
    if (ia < esc_scgrid->na/2)
        return ia + esc_scgrid->na/2;
    else
        return ia - esc_scgrid->na/2;
}

/* Adjust azimuth according to its periodicity: [0;2PI) for the azimuth */
static int sf_esc_scgrid3_bound_ia (sf_esc_scgrid3 esc_scgrid, int ia) {
    if (ia < 0) {
        ia += esc_scgrid->na;
    } else if (ia >= esc_scgrid->na) {
        ia -= esc_scgrid->na;
    }
    return ia;
}

/* Adjust azimuth and inclination angles according to their periodicity:
   [0;PI) for the inclination, [0;2PI) for the azimuth */
static void sf_esc_scgrid3_bound_iaib (sf_esc_scgrid3 esc_scgrid, int *ia, int *ib) {
    /* Azimuth */
    *ia = sf_esc_scgrid3_bound_ia (esc_scgrid, *ia);
    /* Inclination */
    if (*ib < 0) {
        *ib = -(*ib + 1);
        *ia = sf_esc_scgrid3_flip_ia (esc_scgrid, *ia);
    }
    if (*ib >= esc_scgrid->nb) {
        *ib = 2*esc_scgrid->nb - (*ib) - 1;
        *ia = sf_esc_scgrid3_flip_ia (esc_scgrid, *ia);
    }
}

/* Perform bilinear interpolation of the local escape solution */
static void sf_esc_scgrid3_bilinear_interp (sf_esc_scgrid3 esc_scgrid,
                                            float *z, float *x, float *y,
                                            float *t, float *l, float *b, float *a) {
    int i, ib[4], ia[4];
    float vals[4][ESC3_NUM + 3], v[ESC3_NUM + 3];
    float fb, fa, pz, px, py;

    /* ib, ia */
    ib[0] = sf_cram_scgrid3_ftoi (*b, esc_scgrid->ob, esc_scgrid->db, &fb);
    ia[0] = sf_cram_scgrid3_ftoi (*a, esc_scgrid->oa, esc_scgrid->da, &fa);
    /* ib+1, ia */
    ib[1] = ib[0] + 1;
    ia[1] = ia[0];
    /* ib, ia+1 */
    ib[2] = ib[0];
    ia[2] = ia[0] + 1;
    /* ib+1, ia+1 */
    ib[3] = ib[0] + 1;
    ia[3] = ia[2];

    /* Get escape values in space */
    for (i = 0; i < 4; i++) {
        sf_esc_scgrid3_bound_iaib (esc_scgrid, &ia[i], &ib[i]);
        eval_multi_UBspline_3d_s (&esc_scgrid->scsplines[ia[i]*esc_scgrid->nb + ib[i]],
                                  *y, *x, *z, &vals[i][0]);
    }
    /* Bilinear interpolation */
    for (i = 0; i < (ESC3_NUM + 3); i++) {
        v[i] = vals[0][i]*(1.0 - fb)*(1.0 - fa) +
               vals[1][i]*fb*(1.0 - fa) +
               vals[2][i]*(1.0 - fb)*fa +
               vals[3][i]*fb*fa;
    }
    *z += v[ESC3_Z];
    *x += v[ESC3_X];
    *y += v[ESC3_Y];
    *t += v[ESC3_T];
#ifdef ESC_EQ_WITH_L
    *l += v[ESC3_L];
#endif
    pz = cosf (*b) + v[ESC3_NUM];
    px = cosf (*a) + v[ESC3_NUM + 1];
    py = sinf (*a) + v[ESC3_NUM + 2];
    sf_esc_scgrid3_p_to_ab (esc_scgrid, pz, px, py, a, b);
}

#define SPEPS 1e-3

static void sf_esc_scgrid3_project_point (sf_esc_scgrid3 esc_scgrid,
                                          float *z, float *x, float *y,
                                          float *t, float *l, float *b, float *a) {
    if ((*z <= esc_scgrid->zmin + SPEPS*esc_scgrid->dz) ||
        (*z >= esc_scgrid->zmax - SPEPS*esc_scgrid->dz) ||
        (*x <= esc_scgrid->xmin + SPEPS*esc_scgrid->dx) ||
        (*x >= esc_scgrid->xmax - SPEPS*esc_scgrid->dx) ||
        (*y <= esc_scgrid->ymin + SPEPS*esc_scgrid->dy) ||
        (*y >= esc_scgrid->ymax - SPEPS*esc_scgrid->dy)) {
        /* Do ray tracing, if the point is outside of the supergrid bounds */
        sf_esc_tracer3_compute (esc_scgrid->esc_tracer, *z, *x, *y, *b, *a,
                                *t, *l, esc_scgrid->esc_point, a, b);
        *z = sf_esc_point3_get_esc_var (esc_scgrid->esc_point, ESC3_Z);
        *x = sf_esc_point3_get_esc_var (esc_scgrid->esc_point, ESC3_X);
        *y = sf_esc_point3_get_esc_var (esc_scgrid->esc_point, ESC3_Y);
        *t = sf_esc_point3_get_esc_var (esc_scgrid->esc_point, ESC3_T);
#ifdef ESC_EQ_WITH_L
        *l = sf_esc_point3_get_esc_var (esc_scgrid->esc_point, ESC3_L);
#endif
    } else {
        /* Do interpolation of local escape values across the supercells */
        sf_esc_scgrid3_bilinear_interp (esc_scgrid, z, x, y, t, l, b, a);
    }
}

void sf_esc_scgrid3_compute (sf_esc_scgrid3 esc_scgrid,
                             float *z, float *x, float *y,
                             float *t, float *l, float *b, float *a)
/*< Compute escape values for a point with subsurface coordinates (z, x, y, b, a)
    by stitching local escape solutions in supercells of a phase-space grid >*/
{
    while (sf_esc_tracer3_inside (esc_scgrid->esc_tracer, z, x, y, true)) {
        sf_esc_scgrid3_project_point (esc_scgrid, z, x, y, t, l, b, a);
        esc_scgrid->is++;
    }
    esc_scgrid->ir++;
}

