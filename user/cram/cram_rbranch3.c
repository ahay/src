/* Determine surface-exiting ray branches from escape values in 3-D */
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

#include <rsf.h>

#include "esc_point3.h"

#include "cram_slowness3.h"

#ifndef _cram_rbranch3_h

typedef struct CRAMRBranch3 *sf_cram_rbranch3;
/* abstract data type */
/*^*/

typedef struct {
    int   kmah;
    float ib, ia; /* Inclination and azimuth angles */
    float t, cs, s; /* Exit time, exit cosine, exist surface */
    float j; /* Determinant of d(esc_x/y)/d(b/a) */
    float p[3]; /* Slope and its d/dx, d/dy components */
    float ibxy[2], iaxy[2]; /* db/dx, db/dy; da/dx, da/dy */
/*
    float jxy;
*/
} sf_cram_surface_exit3;
/*^*/

#endif

typedef struct {
    int   ib[3], ia[3], kmah; /* Angles, caustic index */
    float xmin, xmax, ymin, ymax; /* Span on the surface */
    float x[3], y[3], t[3]; /* Exit locations and times */
    float p[3], cs, s; /* Slope, slope components (d/dx, d/dy), exit cosine, exit area */
    float j; /* Determinant of d(esc_x/y)/d(x/y) */
    float ibxy[2], iaxy[2]; /* db/dx, db/dy; da/dx, da/dy */
    /* Determinant of d(ib/ia)/d(x/y) */
/*
    float jxy;
*/
} sf_cram_surface_branch3;

struct CRAMRBranch3 {
    int                         na, nb, n, nn, nnx, nny;
    float                       z0, tmax, dbx, dby;
    float                       xmin, xmax, ymin, ymax;
    float                       armin, armax;
    float                       xbmin, xbmax, ybmin, ybmax;
    sf_cram_slowness3           slowness;
    sf_cram_surface_branch3    *surface_branches;
    sf_cram_surface_branch3   **surface_bins;
    sf_cram_surface_branch3  ***surface_ibins;
    int                       **ns;
};
/* concrete data type */

sf_cram_rbranch3 sf_cram_rbranch3_init (int nb, int na, float z0, float tmax,
                                        float xbmin, float xbmax, float ybmin, float ybmax,
                                        float dbx, float dby, sf_cram_slowness3 slowness)
/*< Initialize object >*/
{
    sf_cram_rbranch3 cram_rbranch = (sf_cram_rbranch3)sf_alloc (1, sizeof(struct CRAMRBranch3));

    cram_rbranch->nb = nb; /* Number of escape inclination angles */
    cram_rbranch->na = na; /* Number of escape azimuth angles */
    cram_rbranch->z0 = z0; /* Escape surface (where data are) */
    cram_rbranch->tmax = tmax; /* Maximum time along a ray */

    cram_rbranch->xmin = SF_HUGE;
    cram_rbranch->xmax = -SF_HUGE;
    cram_rbranch->ymin = SF_HUGE;
    cram_rbranch->ymax = -SF_HUGE;

    /* Min & max x, y values for escape coordinates;
       this affects binning and further triangle search
       across the bins */
    cram_rbranch->xbmin = xbmin;
    cram_rbranch->xbmax = xbmax;
    cram_rbranch->ybmin = ybmin;
    cram_rbranch->ybmax = ybmax;
    cram_rbranch->dbx = dbx; /* Surface search bin size in x */
    cram_rbranch->dby = dby; /* Surface search bin size in y */

    cram_rbranch->armax = SF_HUGE;
    cram_rbranch->armin = 1e-5;

    cram_rbranch->surface_branches = NULL;

    /* Number of search bins */
    cram_rbranch->nnx = (int)((cram_rbranch->xbmax - cram_rbranch->xbmin)/
                              cram_rbranch->dbx) + 1;
    cram_rbranch->nny = (int)((cram_rbranch->ybmax - cram_rbranch->ybmin)/
                              cram_rbranch->dby) + 1;
    /* Pointers to the array of search bins */
    cram_rbranch->surface_ibins = (sf_cram_surface_branch3***)sf_alloc (cram_rbranch->nnx*cram_rbranch->nny,
                                                                        sizeof(sf_cram_surface_branch3**));
    /* Array for the number of triangles in each search bin */
    cram_rbranch->ns = sf_intalloc2 (cram_rbranch->nnx, cram_rbranch->nny);
    memset (cram_rbranch->ns[0], 0, cram_rbranch->nnx*cram_rbranch->nny*sizeof(int));
    cram_rbranch->surface_bins = NULL;
    cram_rbranch->nn = 0;

    cram_rbranch->slowness = slowness;

    return cram_rbranch;
}

void sf_cram_rbranch3_close (sf_cram_rbranch3 cram_rbranch)
/*< Destroy object >*/
{
    if (cram_rbranch->surface_bins)
        free (cram_rbranch->surface_bins);
    free (cram_rbranch->surface_ibins);
    free (cram_rbranch->ns[0]);
    free (cram_rbranch->ns);
    if (cram_rbranch)
        free (cram_rbranch->surface_branches);
    free (cram_rbranch);
}

void sf_cram_rbranch3_set_arminmax (sf_cram_rbranch3 cram_rbranch,
                                    float armin, float armax)
/*< Set min and max allowed area for an exit ray on the surface >*/
{
    cram_rbranch->armin = armin;
    cram_rbranch->armax = armax;
}

/*
   Possible triangles (as seen from a,b-coordinate system):
   1----3
   |   /
   |  / 2
   | / /|
   |/ / |
   2 /  |
    /   |
   3----1
   Vectors: 1 = 2-1, 2 = 3-1, 3=2-3
   Vector 1 is in constant azimuth plane
*/

/* Calculate triangle area; db/dx, db/dy and da/dx, da/dy vectors
   around point 1; slope(p) vector around point 1;
   cross product for vectors 1 and 2; and dot product between
   vector 1 and azimuthal vector (defined by angle a) */
static float sf_cram_rbranch3_tridata (float *esc1, float *esc2, float *esc3,
                                       float *ibs, float *ias, float a,
                                       float *ibxy, float *iaxy, float *pxy,
                                       float *cp, float *dp) {
    double s, l1, l2, l3, b1, b2, det;
    double dx1, dy1, dx2, dy2, dx3, dy3; 

    /* Vector 1 */
    dx1 = esc2[ESC3_X] - esc1[ESC3_X];
    dy1 = esc2[ESC3_Y] - esc1[ESC3_Y];
    /* Vector 2 */
    dx2 = esc3[ESC3_X] - esc1[ESC3_X];
    dy2 = esc3[ESC3_Y] - esc1[ESC3_Y];
    /* Vector 3 */
    dx3 = esc2[ESC3_X] - esc3[ESC3_X];
    dy3 = esc2[ESC3_Y] - esc3[ESC3_Y];

    /* Triangle sides length */
    l1 = hypotf (dx1, dy1);
    l2 = hypotf (dx2, dy2);
    l3 = hypotf (dx3, dy3);

    /* 2x2 matrix - normalized components of vectors 1 and 2 */
    det = (dx1*dy2 - dy1*dx2)/(l1*l2);
    if (fabs (det) < 1e-6)
        return -1.0;

    /* Restore p(slope) vector from directional
       derivatives along the sides */
    /* Right-hand side - dt/dl */
    b1 = (esc2[ESC3_T] - esc1[ESC3_T])/l1;
    b2 = (esc3[ESC3_T] - esc1[ESC3_T])/l2;
    /* Slope vector components */
    pxy[0] = (dy2/l2*b1 - dy1/l1*b2)/det;
    pxy[1] = (-dx2/l2*b1 + dx1/l1*b2)/det;

    /* Restore db/dx,db/dy vector from directional
       derivatives along the sides */
    /* Right-hand side - db/dl */
    b1 = (ibs[1] - ibs[0])/l1;
    b2 = (ibs[2] - ibs[0])/l2;
    ibxy[0] = (dy2/l2*b1 - dy1/l1*b2)/det;
    ibxy[1] = (-dx2/l2*b1 + dx1/l1*b2)/det;

    /* Restore da/dx,da/dy vector from directional
       derivatives along the sides */
    /* Right-hand side - da/dl */
    b1 = (ias[1] - ias[0])/l1;
    b2 = (ias[2] - ias[0])/l2;
    iaxy[0] = (dy2/l2*b1 - dy1/l1*b2)/det;
    iaxy[1] = (-dx2/l2*b1 + dx1/l1*b2)/det;

    /* Cross product between side 1 and side 2 */
    *cp = det;

    /* Dot product between side 1 and the azimuthal direction */
    *dp = (cosf (a)*dx1 + sinf (a)*dy1)/l1;

    /* Semi-perimeter */
    s = 0.5*(l1 + l2 + l3);

    /* Return triangle area */
    return sqrt (s*(s - l1)*(s - l2)*(s - l3));
}

static void sf_cram_rbranch3_minmax (float f, float *fmin, float *fmax) {
    if (f < *fmin)
        *fmin = f;
    if (f > *fmax)
        *fmax = f;
}

#define XEPS 1e-4 

static bool sf_cram_rbranch3_check_exit (sf_cram_rbranch3 cram_rbranch,
                                         float *ibs, float *ias, bool left,
                                         float *esc1, float *esc2, float *esc3,
                                         sf_cram_surface_branch3 *branch) {
    int i;
    int kmah = 0;
    float s, a, da, p, sn, vsurf;
    float pxy[2] = {0.0f,0.0f}, cp = 1.0, dp = 0.0;

    /* Check if all three points exit on the surface */
    if (esc1[ESC3_Z] > cram_rbranch->z0 || esc2[ESC3_Z] > cram_rbranch->z0 ||
        esc3[ESC3_Z] > cram_rbranch->z0)
        return false;
    /* Check for min/max time */
    if (esc1[ESC3_T] < 0.0 || esc2[ESC3_T] < 0.0 || esc3[ESC3_T] < 0.0)
        return false;
    if (esc1[ESC3_T] > cram_rbranch->tmax || esc2[ESC3_T] > cram_rbranch->tmax ||
        esc3[ESC3_T] > cram_rbranch->tmax)
        return false;

    da = 2.0*SF_PI/(float)cram_rbranch->na;
    a = 0.5*da + ias[0]*da;
    s = sf_cram_rbranch3_tridata (esc1, esc2, esc3, ibs, ias, a,
                                  branch->ibxy, branch->iaxy,
                                  pxy, &cp, &dp);

    if (s < cram_rbranch->armin || s > cram_rbranch->armax) /* Area is too small or too big */
        return false;

    /* Change sign of the dot product for a left-type triangle */
    if (left)
        dp *= -1.0;

    if (cp >= 0.0) {
        if (fabs (dp) <= 1.0) {
            a = acosf (dp);
            /* Second order caustic - angle between vector
               of constant azimuth and azimuth itself is too high */ 
            if (a > SF_PI/2.0)
                kmah = 2;
        }
    } else /* Negative cross product - first order caustic */
        kmah = 1;

    vsurf = 1.0/sf_cram_slowness3_get_surf_value (cram_rbranch->slowness,
                                                  esc1[ESC3_X], esc1[ESC3_Y]);
    p = hypotf (pxy[0], pxy[1]);
    /* Sine of exit angle */
    sn = fabsf (vsurf*p);
    if (sn >= 1.0) {
        /* Too wide ray tube */
        return false;
    }
    branch->s = s;
    /* Exit cosine */
    branch->cs = sqrt (1.0 - sn*sn);
    /* Determinant of d(esc_x/y)/d(b/a) */
    branch->j = (esc2[ESC3_X] - esc1[ESC3_X])/(float)(ibs[1] - ibs[0])*
                (esc3[ESC3_Y] - esc1[ESC3_Y])/(float)(ias[2] - ias[0]) -
                (esc2[ESC3_Y] - esc1[ESC3_Y])/(float)(ibs[1] - ibs[0])*
                (esc3[ESC3_X] - esc1[ESC3_X])/(float)(ias[2] - ias[0]);
    /* Determinant of d(ib/ia)/d(x/y) */
/*
    branch->jxy = branch->ibxy[0]*branch->iaxy[1] - 
                  branch->ibxy[1]*branch->iaxy[0];
*/
    branch->kmah = kmah;
    for (i = 0; i < 3; i++) { 
        branch->ib[i] = ibs[i];
        branch->ia[i] = ias[i];
    }
    branch->x[0] = esc1[ESC3_X]; branch->y[0] = esc1[ESC3_Y]; branch->t[0] = esc1[ESC3_T];
    branch->x[1] = esc2[ESC3_X]; branch->y[1] = esc2[ESC3_Y]; branch->t[1] = esc2[ESC3_T];
    branch->x[2] = esc3[ESC3_X]; branch->y[2] = esc3[ESC3_Y]; branch->t[2] = esc3[ESC3_T];
    branch->p[0] = p;
    branch->p[1] = pxy[0];
    branch->p[2] = pxy[1];

    branch->xmin = SF_HUGE;
    branch->xmax = -SF_HUGE;
    branch->ymin = SF_HUGE;
    branch->ymax = -SF_HUGE;

    /* Local min/max */
    sf_cram_rbranch3_minmax (esc1[ESC3_X], &branch->xmin, &branch->xmax);
    sf_cram_rbranch3_minmax (esc1[ESC3_Y], &branch->ymin, &branch->ymax);
    sf_cram_rbranch3_minmax (esc2[ESC3_X], &branch->xmin, &branch->xmax);
    sf_cram_rbranch3_minmax (esc2[ESC3_Y], &branch->ymin, &branch->ymax);
    sf_cram_rbranch3_minmax (esc3[ESC3_X], &branch->xmin, &branch->xmax);
    sf_cram_rbranch3_minmax (esc3[ESC3_Y], &branch->ymin, &branch->ymax);
    /* Global min/max */
    sf_cram_rbranch3_minmax (esc1[ESC3_X], &cram_rbranch->xmin, &cram_rbranch->xmax);
    sf_cram_rbranch3_minmax (esc1[ESC3_Y], &cram_rbranch->ymin, &cram_rbranch->ymax);
    sf_cram_rbranch3_minmax (esc2[ESC3_X], &cram_rbranch->xmin, &cram_rbranch->xmax);
    sf_cram_rbranch3_minmax (esc2[ESC3_Y], &cram_rbranch->ymin, &cram_rbranch->ymax);
    sf_cram_rbranch3_minmax (esc3[ESC3_X], &cram_rbranch->xmin, &cram_rbranch->xmax);
    sf_cram_rbranch3_minmax (esc3[ESC3_Y], &cram_rbranch->ymin, &cram_rbranch->ymax);

    return true;
}

void sf_cram_rbranch3_set_escapes (sf_cram_rbranch3 cram_rbranch, float ***esc)
/*< Set new escape variables >*/
{
    int ib, ia, i, iam, iix, iiy;
    int iix0, iix1, iiy0, iiy1;
    float ibs[3], ias[3], cap[ESC3_NUM];
    sf_cram_surface_branch3 *branch, *branches;
    sf_cram_surface_branch3 **bin;

    if (cram_rbranch->surface_branches)
        free (cram_rbranch->surface_branches);
    cram_rbranch->surface_branches = (sf_cram_surface_branch3*)sf_alloc (cram_rbranch->na*cram_rbranch->nb*2,
                                                                         sizeof(sf_cram_surface_branch3));

    cram_rbranch->xmin = SF_HUGE;
    cram_rbranch->xmax = -SF_HUGE;
    cram_rbranch->ymin = SF_HUGE;
    cram_rbranch->ymax = -SF_HUGE;
    cram_rbranch->n = 0;

    /* Approximate escape value at the pole */
    cap[ESC3_Z] = 0.5*(esc[0][0][ESC3_Z] + esc[cram_rbranch->na/2][0][ESC3_Z]);
    cap[ESC3_X] = 0.5*(esc[0][0][ESC3_X] + esc[cram_rbranch->na/2][0][ESC3_X]);
    cap[ESC3_Y] = 0.5*(esc[0][0][ESC3_Y] + esc[cram_rbranch->na/2][0][ESC3_Y]);
    cap[ESC3_T] = 0.5*(esc[0][0][ESC3_T] + esc[cram_rbranch->na/2][0][ESC3_T]);

    /* Scan ray branches and store those exiting on the surface */
    for (ia = 0; ia < cram_rbranch->na; ia++) {
        iam = ia != 0 ? ia - 1 : cram_rbranch->na - 1;
        /* Process the pole cap */
        ibs[0] = 0; ibs[1] = -0.5; ibs[2] = 0;
        ias[0] = ia; ias[1] = ia - 0.5; ias[2] = ia - 1;
        cram_rbranch->n += sf_cram_rbranch3_check_exit (cram_rbranch, ibs, ias, true,
                                                        esc[ia][0], cap, esc[iam][0],
                                                        &cram_rbranch->surface_branches[cram_rbranch->n]);
        /* Process remaining triangles in the azimuth band */
        for (ib = 1; ib < cram_rbranch->nb; ib++) {
            /* Break every escape polygon into a pair of triangles */
            ibs[0] = ib - 1; ibs[1] = ib; ibs[2] = ib - 1;
            ias[0] = ia - 1; ias[1] = ia - 1; ias[2] = ia;
            cram_rbranch->n += sf_cram_rbranch3_check_exit (cram_rbranch, ibs, ias, false,
                                                            esc[iam][ib - 1], esc[iam][ib], esc[ia][ib - 1],
                                                            &cram_rbranch->surface_branches[cram_rbranch->n]);
            ibs[0] = ib; ibs[1] = ib - 1; ibs[2] = ib;
            ias[0] = ia; ias[1] = ia; ias[2] = ia - 1;
            cram_rbranch->n += sf_cram_rbranch3_check_exit (cram_rbranch, ibs, ias, true,
                                                            esc[ia][ib], esc[ia][ib - 1], esc[iam][ib],
                                                            &cram_rbranch->surface_branches[cram_rbranch->n]);
        }
    }
    memset (cram_rbranch->ns[0], 0, cram_rbranch->nnx*cram_rbranch->nny*sizeof(int));
    cram_rbranch->nn = 0;
    if (0 == cram_rbranch->n)
        return;

    /* Reallocate branches storage to the actual useful size */
    branches = (sf_cram_surface_branch3*)sf_alloc (cram_rbranch->n, sizeof(sf_cram_surface_branch3));
    memcpy (branches, cram_rbranch->surface_branches, sizeof(sf_cram_surface_branch3)*cram_rbranch->n);
    free (cram_rbranch->surface_branches);
    cram_rbranch->surface_branches = branches;

    /* Find out how many triangles are in each search bin */
    for (i = 0; i < cram_rbranch->n; i++) {
        branch = &cram_rbranch->surface_branches[i];
        iix0 = (int)((branch->xmin - cram_rbranch->xbmin)/cram_rbranch->dbx);
        iix1 = (int)((branch->xmax - cram_rbranch->xbmin)/cram_rbranch->dbx);
        if (iix0 < 0)
            iix0 = 0;
        if (iix1 >= cram_rbranch->nnx)
            iix1 = cram_rbranch->nnx - 1;
        iiy0 = (int)((branch->ymin - cram_rbranch->ybmin)/cram_rbranch->dby);
        iiy1 = (int)((branch->ymax - cram_rbranch->ybmin)/cram_rbranch->dby);
        if (iiy0 < 0)
            iiy0 = 0;
        if (iiy1 >= cram_rbranch->nny)
            iiy1 = cram_rbranch->nny - 1;
        for (iiy = iiy0; iiy <= iiy1; iiy++) {
            for (iix = iix0; iix <= iix1; iix++) {
                cram_rbranch->ns[iiy][iix]++;
                cram_rbranch->nn++;
            }
        }
    }

    if (cram_rbranch->surface_bins)
        free (cram_rbranch->surface_bins);
    cram_rbranch->surface_bins = (sf_cram_surface_branch3**)sf_alloc (cram_rbranch->nn,
                                                                      sizeof(sf_cram_surface_branch3*));
    /* Create indices for search bins */
    i = 0;
    for (iiy = 0; iiy < cram_rbranch->nny; iiy++) {
        for (iix = 0; iix < cram_rbranch->nnx; iix++) {
            cram_rbranch->surface_ibins[iiy*cram_rbranch->nnx + iix] = &cram_rbranch->surface_bins[i];
            i += cram_rbranch->ns[iiy][iix];
        }
    }

    /* Fill the bins with pointers to the triangles */
    memset (cram_rbranch->ns[0], 0, cram_rbranch->nnx*cram_rbranch->nny*sizeof(int));
    for (i = 0; i < cram_rbranch->n; i++) {
        branch = &cram_rbranch->surface_branches[i];
        iix0 = (int)((branch->xmin - cram_rbranch->xbmin)/cram_rbranch->dbx);
        iix1 = (int)((branch->xmax - cram_rbranch->xbmin)/cram_rbranch->dbx);
        if (iix0 < 0)
            iix0 = 0;
        if (iix1 >= cram_rbranch->nnx)
            iix1 = cram_rbranch->nnx - 1;
        iiy0 = (int)((branch->ymin - cram_rbranch->ybmin)/cram_rbranch->dby);
        iiy1 = (int)((branch->ymax - cram_rbranch->ybmin)/cram_rbranch->dby);
        if (iiy0 < 0)
            iiy0 = 0;
        if (iiy1 >= cram_rbranch->nny)
            iiy1 = cram_rbranch->nny - 1;
        for (iiy = iiy0; iiy <= iiy1; iiy++) {
            for (iix = iix0; iix <= iix1; iix++) {
                bin = cram_rbranch->surface_ibins[iiy*cram_rbranch->nnx + iix];
                bin[cram_rbranch->ns[iiy][iix]] = branch;
                cram_rbranch->ns[iiy][iix]++;
            }
        }
    }
}

/* Returns barycentric coordinates (u,v) for a point (px,py) in a triangle abc */
static bool sf_cram_rbranch3_barycentric (float ax, float ay, float bx, float by,
                                          float cx, float cy, float px, float py,
                                          double *u, double *v) {
   double v0x = cx - ax, v0y = cy - ay,
         v1x = bx - ax, v1y = by - ay,
         v2x = px - ax, v2y = py - ay;
   double dot00 = v0x*v0x + v0y*v0y,
         dot01 = v0x*v1x + v0y*v1y,
         dot02 = v0x*v2x + v0y*v2y,
         dot11 = v1x*v1x + v1y*v1y,
         dot12 = v1x*v2x + v1y*v2y;
   double invDenom = 1.0 / (dot00*dot11 - dot01*dot01);

   *u = (dot11*dot02 - dot01*dot12)*invDenom;
   *v = (dot00*dot12 - dot01*dot02)*invDenom;

   return (*u >= 0.0) && (*v >= 0.0) && (*u + *v <= 1.0);
}

int sf_cram_rbranch3_find_exits (sf_cram_rbranch3 cram_rbranch, float x, float y,
                                 int ne, sf_cram_surface_exit3* exits)
/*< Find all exit triangles for the point x, y on the surface
    and fill exits[ne] array with pointers to them >*/
{
    int i, ie, iie, iix, iiy, n;
    float t;
    bool same;
    double u, v;
    sf_cram_surface_branch3 *branch;
    sf_cram_surface_branch3 **bin;

    if (x > cram_rbranch->xmax || x < cram_rbranch->xmin ||
        y > cram_rbranch->ymax || y < cram_rbranch->ymin)
        return 0;

    i = 0;
    ie = 0;

    /* Isolate search bin */
    iix = (int)((x - cram_rbranch->xbmin)/cram_rbranch->dbx);
    iiy = (int)((y - cram_rbranch->ybmin)/cram_rbranch->dby);
    bin = cram_rbranch->surface_ibins[iiy*cram_rbranch->nnx + iix];
    n = cram_rbranch->ns[iiy][iix];
    /* Check for triangle hits within the search bin */
    while (i < n && ie < ne) {
        branch = bin[i];
        if (x >= branch->xmin && x <= branch->xmax &&
            y >= branch->ymin && y <= branch->ymax &&
            sf_cram_rbranch3_barycentric (branch->x[0], branch->y[0],
                                          branch->x[1], branch->y[1],
                                          branch->x[2], branch->y[2],
                                          x, y, &u, &v)) {
            /* Linear interpolation of traveltime from 3 triangle vertices */ 
            t = branch->t[0]*(1.0 - u - v) + branch->t[1]*v + branch->t[2]*u;
            if (t <= cram_rbranch->tmax) {
                exits[ie].ib = branch->ib[0]*(1.0 - u - v) + branch->ib[1]*v + branch->ib[2]*u;
                exits[ie].ia = branch->ia[0]*(1.0 - u - v) + branch->ia[1]*v + branch->ia[2]*u;
                same = false;
                /* Check if this arrival is different from those previously found;
                   duplication of arrivals can happen if very small exit triangles
                   are not filtered out */
                for (iie = 0; iie < ie; iie++) {
                    if (fabsf (t - exits[iie].t) < 1e-3 && /* Same traveltime */
                        exits[ie].kmah == branch->kmah &&  /* and same caustic */
                        ((fabsf (exits[ie].ib - exits[iie].ib) < 1e-3 && /* + Same angles */
                          fabsf (exits[ie].ia - exits[iie].ia) < 1e-3) ||
                         fabsf (branch->p[0] - exits[iie].p[0]) < 1e-2 || /* or same data slope */
                         (exits[ie].ib < 0.5 && exits[iie].ib < 0.5))) {  /* or close to the pole */
                        same = true;
                        break;
                    }
                }
                if (!same) {
                    exits[ie].t = t;
                    exits[ie].j = branch->j;
/*
                    exits[ie].jxy = branch->jxy;
*/
                    exits[ie].s = branch->s;
                    exits[ie].cs = branch->cs;
                    exits[ie].p[0] = branch->p[0];
                    exits[ie].p[1] = branch->p[1];
                    exits[ie].p[2] = branch->p[2];
                    exits[ie].ibxy[0] = branch->ibxy[0];
                    exits[ie].ibxy[1] = branch->ibxy[1];
                    exits[ie].iaxy[0] = branch->iaxy[0];
                    exits[ie].iaxy[1] = branch->iaxy[1];
                    exits[ie].kmah = branch->kmah;
                    ie++;
                }
            }
        }
        i++;
    }

    return ie;
}

