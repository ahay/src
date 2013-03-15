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

typedef struct {
    int   iab;
    float z, x, y;
} sf_esc_scgrid3_areq;
/* Structure for requesting one (z,x,y) point in angle space */
/*^*/

#endif

#include "einspline.h"
#include "esc_helper.h"

struct EscSCgrid3 {
    size_t               n, offs, is, ir;
    int                  morder;
    int                  nz, nx, ny, na, nb;
    float                oz, ox, oy, oa, ob;
    float                dz, dx, dy, da, db;
    float                zmin, zmax, xmin, xmax, ymin, ymax, md;
    bool                 quad;
    multi_UBspline_3d_s *scsplines;
    unsigned char       *mmaped;
    sf_esc_point3        esc_point;
    sf_esc_tracer3       esc_tracer;
    float              **L;
    int                 *pvt, ns;
    float               *bs, *as;
};
/* concrete data type */

/* 8-point stencil, version one:
   b
   ^
   |
   +---+---X---+
   |   |   |   |
   X---X---X---+
   |   |   |   |
   +---X---x---X
   |   |   |   |
   +---X---+---+--->a
*/
/*
#define SCGRID3_TPS_STENCIL1 8
static int sf_scgrid3_tps_ib_stencil1[SCGRID3_TPS_STENCIL1] = 
{ 2, 0, 1, 2, 1, 2, 3, 1 };
static int sf_scgrid3_tps_ia_stencil1[SCGRID3_TPS_STENCIL1] = 
{ 0, 1, 1, 1, 2, 2, 2, 3 };
*/
/* 8-point stencil, version two:
   b
   ^
   |
   X---+---+---X
   |   |   |   |
   +---X---X---+
   |   |   |   |
   +---X---x---+
   |   |   |   |
   X---+---+---X--->a
*/
/*
#define SCGRID3_TPS_STENCIL2 8
static int sf_scgrid3_tps_ib_stencil2[SCGRID3_TPS_STENCIL2] = 
{ 0, 3, 1, 2, 1, 2, 0, 3 };
static int sf_scgrid3_tps_ia_stencil2[SCGRID3_TPS_STENCIL2] = 
{ 0, 0, 1, 1, 2, 2, 3, 3 };
*/
/* 12-point stencil:
   b
   ^
   |
   +---X---X---+
   |   |   |   |
   X---X---X---X
   |   |   |   |
   X---X---x---X
   |   |   |   |
   +---X---X---+--->a
*/

#define SCGRID3_TPS_STENCIL3 12
static int sf_scgrid3_tps_ib_stencil3[SCGRID3_TPS_STENCIL3] = 
{ 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 1, 2 };
static int sf_scgrid3_tps_ia_stencil3[SCGRID3_TPS_STENCIL3] = 
{ 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3 };

/* 16-point stencil:
   b
   ^
   |
   X---X---X---X
   |   |   |   |
   X---X---X---X
   |   |   |   |
   X---X---x---X
   |   |   |   |
   X---X---X---X--->a
*/
/*
#define SCGRID3_TPS_STENCIL4 16
static int sf_scgrid3_tps_ib_stencil4[SCGRID3_TPS_STENCIL4] = 
{ 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3 };
static int sf_scgrid3_tps_ia_stencil4[SCGRID3_TPS_STENCIL4] = 
{ 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3 };
*/
#define SCGRID3_TPS_MAX_STENCIL 16

/* Initialize thin-plane spline structures */
static void sf_esc_scgrid3_init_tps (sf_esc_scgrid3 esc_scgrid) {
    const int ns = SCGRID3_TPS_STENCIL3;
    int *bs = sf_scgrid3_tps_ib_stencil3;
    int *as = sf_scgrid3_tps_ia_stencil3;
    int i;

    esc_scgrid->L = sf_floatalloc2 (ns + 3, ns + 3);
    esc_scgrid->as = sf_floatalloc (ns);
    esc_scgrid->bs = sf_floatalloc (ns);
    esc_scgrid->pvt = sf_intalloc (ns + 3);

    for (i = 0; i < ns; i++) {
        esc_scgrid->bs[i] = (float)bs[i];
        esc_scgrid->as[i] = (float)as[i];
    }
    sf_tps_build_matrix (esc_scgrid->L, esc_scgrid->bs, esc_scgrid->as, ns, 0.0);
    if (0 != sf_ludlt_decomposition (&esc_scgrid->L[0][0], esc_scgrid->pvt,            
                                     ns + 3)) {
        sf_error ("sf_esc_scgrid3_init_tps: TPS matrix inverse failed"); 
    }
    esc_scgrid->ns = ns;
}

sf_esc_scgrid3 sf_esc_scgrid3_init (sf_file scgrid, sf_esc_tracer3 esc_tracer, bool verb)
/*< Initialize object >*/
{
    size_t nc, nnc = 0;
    int ia, ib;
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
    if (!sf_histint (scgrid, "Nb", &esc_scgrid->nb)) sf_error ("No Nb= in supercell file");
    if (!sf_histfloat (scgrid, "Db", &esc_scgrid->db)) sf_error ("No Db= in supercell file");
    if (!sf_histfloat (scgrid, "Ob", &esc_scgrid->ob)) sf_error ("No Ob= in supercell file");
    if (!sf_histint (scgrid, "Na", &esc_scgrid->na)) sf_error ("No Na= in supercell file");
    if (!sf_histfloat (scgrid, "Da", &esc_scgrid->da)) sf_error ("No Da= in supercell file");
    if (!sf_histfloat (scgrid, "Oa", &esc_scgrid->oa)) sf_error ("No Oa= in supercell file");

    if (!sf_histfloat (scgrid, "Mdist", &esc_scgrid->md)) esc_scgrid->md = esc_scgrid->dz;

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
#ifdef DEBUG
    esc_scgrid->mmaped = sf_ucharalloc ((size_t)esc_scgrid->offs +
                                        (size_t)esc_scgrid->n*
                                        (size_t)esc_scgrid->na*
                                        (size_t)esc_scgrid->nb);
    sf_ucharread (esc_scgrid->mmaped, (size_t)esc_scgrid->offs +
                                      (size_t)esc_scgrid->n*
                                      (size_t)esc_scgrid->na*
                                      (size_t)esc_scgrid->nb, scgrid);
#else
    esc_scgrid->mmaped = (unsigned char*)mmap (NULL, (size_t)esc_scgrid->offs +
                                                     (size_t)esc_scgrid->n*
                                                     (size_t)esc_scgrid->na*
                                                     (size_t)esc_scgrid->nb,
                                               PROT_READ, MAP_SHARED,
                                               fileno (stream), 0);
#endif
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

    esc_scgrid->morder = 2; /* Interpolation accuracy */

    esc_scgrid->esc_tracer = esc_tracer;
    sf_esc_tracer3_set_mdist (esc_tracer, esc_scgrid->md);
    esc_scgrid->esc_point = sf_esc_point3_init ();

    sf_esc_scgrid3_init_tps (esc_scgrid);

    return esc_scgrid;
}

void sf_esc_scgrid3_close (sf_esc_scgrid3 esc_scgrid, bool verb)
/*< Destroy object >*/
{
    if (verb)
        sf_warning ("%lu points processed, %g interpolation steps per point performed",
                    esc_scgrid->ir, (float)esc_scgrid->is/(float)esc_scgrid->ir);
    free (esc_scgrid->L[0]);
    free (esc_scgrid->L);
    free (esc_scgrid->as);
    free (esc_scgrid->bs);
    free (esc_scgrid->pvt);
#ifdef DEBUG
    free (esc_scgrid->mmaped);
#else
    munmap (esc_scgrid->mmaped, (size_t)esc_scgrid->offs +
                                (size_t)esc_scgrid->n*
                                (size_t)esc_scgrid->na*
                                (size_t)esc_scgrid->nb);
#endif
    free (esc_scgrid->scsplines);
    sf_esc_point3_close (esc_scgrid->esc_point);
    free (esc_scgrid);
}

void sf_esc_scgrid3_set_morder (sf_esc_scgrid3 esc_scgrid, int morder)
/*< Set order of interpolation accuracy in the angular domain >*/
{
    esc_scgrid->morder = morder;
}

/* Compute interger index for a float value v according to sampling df
   and zero at f0, store fractional part in f */
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
static void sf_esc_scgrid3_p_to_ab (float pz, float px, float py,
                                    float *a, float *b) {
    float l;

    l = sqrt (px*px + py*py);
    if (l > 1e-6) {
        if (py >= 0.0)
            *a = acosf (px/l);
        else
            *a = 2.0*SF_PI - acosf (px/l);
    }
    l = sqrt (pz*pz + px*px + py*py);
    *b = acosf (pz/l);
}

/* Convert quaternion of rotation into a new phase direction */
static void sf_esc_scgrid3_q_to_ab (float *q, float *a, float *b) {
    int i, j;
    float vf[3], vt[3];
    float M[3][3];

    vf[0] = cosf (*b);
    vf[1] = sinf (*b)*cosf (*a);
    vf[2] = sinf (*b)*sinf (*a);

    /* Obtain rotation to the new phase position */
    sf_quat_norm (q);
    sf_quat_rotmat (q, &M[0][0]);

    /* Rotate the old phase vector */
    for (i = 0; i < 3; i++) {
        vt[i] = 0.0;
        for (j = 0; j < 3; j++) {
            vt[i] += M[i][j]*vf[j];
        }
    }
    /* Convert to azimuth and inclination */
    sf_esc_scgrid3_p_to_ab (vt[0], vt[1], vt[2], a, b);
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

/* Perform nearest neighbor interpolation of the local escape solution */
static void sf_esc_scgrid3_nearest_interp (sf_esc_scgrid3 esc_scgrid,
                                           float *z, float *x, float *y,
                                           float *t, float *l, float *b, float *a) {
    int ib, ia;
    float vals[ESC3_NUM + 4];
    float fb, fa;

    ib = sf_cram_scgrid3_ftoi (*b, esc_scgrid->ob, esc_scgrid->db, &fb);
    ia = sf_cram_scgrid3_ftoi (*a, esc_scgrid->oa, esc_scgrid->da, &fa);

    ib = (int)((float)ib + fb + 0.5);
    ia = (int)((float)ia + fa + 0.5);

    /* Get escape values in space */
    sf_esc_scgrid3_bound_iaib (esc_scgrid, &ia, &ib);
    eval_multi_UBspline_3d_s (&esc_scgrid->scsplines[ia*esc_scgrid->nb + ib],
                              *y, *x, *z, vals);
    *z += vals[ESC3_Z];
    *x += vals[ESC3_X];
    *y += vals[ESC3_Y];
    *t += vals[ESC3_T];
#ifdef ESC_EQ_WITH_L
    *l += vals[ESC3_L];
#endif
    sf_esc_scgrid3_q_to_ab (&vals[ESC3_NUM], a, b);
}

/* Perform bilinear interpolation of the local escape solution */
static void sf_esc_scgrid3_bilinear_interp (sf_esc_scgrid3 esc_scgrid,
                                            float *z, float *x, float *y,
                                            float *t, float *l, float *b, float *a) {
    int i, ib[4], ia[4];
    float vals[4][ESC3_NUM + 4], v[ESC3_NUM + 4];
    float fb, fa;

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
    for (i = 0; i < (ESC3_NUM + 4); i++) {
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
    sf_esc_scgrid3_q_to_ab (&v[ESC3_NUM], a, b);
}

/* Perform spline interpolation of the local escape solution */
static void sf_esc_scgrid3_spline_interp (sf_esc_scgrid3 esc_scgrid,
                                          float *z, float *x, float *y,
                                          float *t, float *l, float *b, float *a) {
    int i, j, k, ia, ib, iia, iib;
    float vals[6][6][ESC3_NUM + 4], f[ESC3_NUM + 4][6][6], v[ESC3_NUM + 3];
    float fa, fb;
    Ugrid z_grid, x_grid;
    BCtype_s zBC, xBC;
    multi_UBspline_2d_s *escspline = NULL;

    /* ib, ia */
    iib = sf_cram_scgrid3_ftoi (*b, esc_scgrid->ob, esc_scgrid->db, &fb);
    iia = sf_cram_scgrid3_ftoi (*a, esc_scgrid->oa, esc_scgrid->da, &fa);

    /* Collect 6x6 array of escape values for 2-D spline interpolation */
    for (i = 0; i < 6; i++) {
        for (j = 0; j < 6; j++) {
            ia = iia + i - 2;
            ib = iib + j - 2;
            sf_esc_scgrid3_bound_iaib (esc_scgrid, &ia, &ib);
            eval_multi_UBspline_3d_s (&esc_scgrid->scsplines[ia*esc_scgrid->nb + ib],
                                      *y, *x, *z, &vals[i][j][0]);
        }
    }
    for (k = 0; k < (ESC3_NUM + 4); k++) {
        for (i = 0; i < 6; i++) {
            for (j = 0; j < 6; j++) {
                f[k][i][j] = vals[i][j][k];
            }
        }
    }
    /* Spline interpolation */
    z_grid.start = 0; z_grid.end = 5; z_grid.num = 6;
    x_grid.start = 0; x_grid.end = 5; x_grid.num = 6;
    zBC.lCode = zBC.rCode = NATURAL;
    xBC.lCode = xBC.rCode = NATURAL;
    escspline = create_multi_UBspline_2d_s (x_grid, z_grid, xBC, zBC, ESC3_NUM + 4);
    for (k = 0; k < (ESC3_NUM + 4); k++) {
        set_multi_UBspline_2d_s (escspline, k, &f[k][0][0]);
    }
    eval_multi_UBspline_2d_s (escspline, 2.0 + fa, 2.0 + fb, v);
    *z += v[ESC3_Z];
    *x += v[ESC3_X];
    *y += v[ESC3_Y];
    *t += v[ESC3_T];
#ifdef ESC_EQ_WITH_L
    *l += v[ESC3_L];
#endif
    sf_esc_scgrid3_q_to_ab (&v[ESC3_NUM], a, b);
    destroy_Bspline (escspline);
}

/* Perform thin plate spline interpolation of the local escape solution */
static void sf_esc_scgrid3_tps_interp (sf_esc_scgrid3 esc_scgrid,
                                       float *z, float *x, float *y,
                                       float *t, float *l, float *b, float *a) {
    int i, j, ia, ib, iia, iib;
    float vals[SCGRID3_TPS_MAX_STENCIL][ESC3_NUM + 4],
          f[SCGRID3_TPS_MAX_STENCIL + 3], w[SCGRID3_TPS_MAX_STENCIL + 3],
          v[ESC3_NUM + 4];
    float fa, fb;

    /* ib, ia */
    iib = sf_cram_scgrid3_ftoi (*b, esc_scgrid->ob, esc_scgrid->db, &fb);
    iia = sf_cram_scgrid3_ftoi (*a, esc_scgrid->oa, esc_scgrid->da, &fa);

    /* Extract escape values according to the interpolation stencil */
    for (i = 0; i < esc_scgrid->ns; i++) {
        ib = esc_scgrid->bs[i] - 1 + iib;
        ia = esc_scgrid->as[i] - 1 + iia;
        sf_esc_scgrid3_bound_iaib (esc_scgrid, &ia, &ib);
        eval_multi_UBspline_3d_s (&esc_scgrid->scsplines[ia*esc_scgrid->nb + ib],
                                  *y, *x, *z, &vals[i][0]);
        
    }
    /* Loop over escape functions */
    for (i = 0; i < (ESC3_NUM + 4); i++) {
        /* Put one escape function type into a contiguous array */ 
        for (j = 0; j < esc_scgrid->ns; j++) {
            f[j] = vals[j][i];
        }
        for (j = esc_scgrid->ns; j < (esc_scgrid->ns + 3); j++) {
            f[j] = 0.0;
        }
        /* Compute thin-plate spline coefficients */
        if (0 != sf_ludtl_solve (&esc_scgrid->L[0][0], f, esc_scgrid->pvt,
                                 w, esc_scgrid->ns + 3)) {
            sf_error ("sf_esc_scgrid3_tps_interp: linear system solution error");
            exit (-1);
        }
        /* Do interpolation */
        v[i] = sf_tps_compute_point (w, esc_scgrid->bs, esc_scgrid->as,
                                     esc_scgrid->ns, 1.0 + fb, 1.0 + fa);
    }
    *z += v[ESC3_Z];
    *x += v[ESC3_X];
    *y += v[ESC3_Y];
    *t += v[ESC3_T];
#ifdef ESC_EQ_WITH_L
    *l += v[ESC3_L];
#endif
    sf_esc_scgrid3_q_to_ab (&v[ESC3_NUM], a, b);
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
        if (3 == esc_scgrid->morder)
            sf_esc_scgrid3_spline_interp (esc_scgrid, z, x, y, t, l, b, a);
        else if (2 == esc_scgrid->morder)
            sf_esc_scgrid3_tps_interp (esc_scgrid, z, x, y, t, l, b, a);
        else if (0 == esc_scgrid->morder)
            sf_esc_scgrid3_nearest_interp (esc_scgrid, z, x, y, t, l, b, a);
        else
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

