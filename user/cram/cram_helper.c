/* Auxiliary functions for angle migration */
/*
  Copyright (C) 2011 University of Texas at Austin

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

#define DHORDER 3

/* F-D coefficients (FDo9p) from Bogey and Bailly, 2004, JCP */
#define C1 0.841570125482
#define C2 -0.244678631765
#define C3 0.059463584768
#define C4 -0.007650904064

static float sf_cram_bbfd_stencil (float am4, float am3, float am2, float am1,
                                   float ap1, float ap2, float ap3, float ap4) {
    return -C4*am4 -C3*am3 -C2*am2 -C1*am1
           +C1*ap1 +C2*ap2 +C3*ap3 +C4*ap4;
}

void sf_cram_trace_deriv (float *tracein, float *traceout, int n, float d)
/*< First order derivative of a trace[n] >*/
{
    int i;
    float invd = 1.0/d;
    float b = tracein[0], e = tracein[n - 1];
    float bd = tracein[0] - tracein[1];
    float ed = tracein[n - 1] - tracein[n - 2];

    for (i = 0; i < n; i++) {
        if (i < 4)
            traceout[i] = sf_cram_bbfd_stencil (b + (4.0 - (float)i)*bd,
                                                i > 2 ? tracein[i - 3] : b + (3.0 - (float)i)*bd,
                                                i > 1 ? tracein[i - 2] : b + (2.0 - (float)i)*bd,
                                                i > 0 ? tracein[i - 1] : b + (1.0 - (float)i)*bd,
                                                tracein[i + 1], tracein[i + 2],
                                                tracein[i + 3], tracein[i + 4]);
        else if (i > (n - 5))
            traceout[i] = sf_cram_bbfd_stencil (tracein[i - 4], tracein[i - 3],
                                                tracein[i - 2], tracein[i - 1],
                                                i < (n - 1) ? tracein[i + 1] :
                                                e + (1.0 - (float)(n - i - 1))*ed,
                                                i < (n - 2) ? tracein[i + 2] :
                                                e + (2.0 - (float)(n - i - 1))*ed,
                                                i < (n - 3) ? tracein[i + 3] :
                                                e + (3.0 - (float)(n - i - 1))*ed,
                                                e + (4.0 - (float)(n - i - 1))*ed);
        else
            traceout[i] = sf_cram_bbfd_stencil (tracein[i - 4], tracein[i - 3],
                                                tracein[i - 2], tracein[i - 1],
                                                tracein[i + 1], tracein[i + 2],
                                                tracein[i + 3], tracein[i + 4]);
        traceout[i] *= invd;
    }
}

void sf_cram_trace_cint (float *trace, int n)
/*< Causal integration of a trace[n] >*/
{
    int i;

    for (i = 1; i < n; i++)
        trace[i] += trace[i - 1];
}

void sf_cram_trace_acint (float *trace, int n)
/*< Anticausal integrations of a trace[n] >*/
{
    int i;

    for (i = n - 2; i >= 0; i--)
        trace[i] += trace[i + 1];
}

/* Convolution operator - borrowed from SU */
static void sf_cram_trace_conv (int lx, int ifx, float *x,
                                int ly, int ify, float *y,
                                int lz, int ifz, float *z) {
    int ilx = ifx + lx - 1,
        ily = ify + ly - 1,
        ilz = ifz + lz - 1,
        i, j, jlow, jhigh;
    float sum;

    x -= ifx;  y -= ify;  z -= ifz;
    for (i=ifz; i<=ilz; ++i) {
        jlow = i - ily;
        if (jlow < ifx)
            jlow = ifx;
        jhigh = i-ify;
        if (jhigh > ilx)
            jhigh = ilx;
        for (j = jlow, sum = 0.0; j <= jhigh; ++j)
            sum += x[j]*y[i-j];
        z[i] = sum;
    }
}

/* Hilbert transform - borrowed from SU */
#define LHHALF 30       /* half-length of Hilbert transform filter*/
#define LH 2*LHHALF+1   /* filter length must be odd */
void sf_cram_trace_hilbert (int n, float *x, float *y)
/*< Hilbert tracnform of a trace x[n] -> y[n] >*/
{
    static int madeh=0;
    static float h[LH];
    int i;
    float taper;

    /* if not made, make Hilbert transform filter; use Hamming window */
    if (!madeh) {
        h[LHHALF] = 0.0;
        for (i=1; i<=LHHALF; i++) {
            taper = 0.54 + 0.46*cos (M_PI*(float)i/(float)(LHHALF));
            h[LHHALF+i] = taper*(-(float)(i % 2)*2.0/
                         (M_PI*(float)(i)));
            h[LHHALF-i] = -h[LHHALF + i];
        }
        madeh = 1;
    }

    /* convolve Hilbert transform with input array */
    sf_cram_trace_conv (LH, -LHHALF, h, n, 0, x, n, 0, y);
}

/*
 *  2-D Delaunay triangulation based of Paul Bourke's code
 *  
 */

#define TRI_EPS 1e-6
 
/*
   Returns true, if a point (xp,yp) is inside the circumcircle made up
   of the points (x1,y1), (x2,y2), (x3,y3).
   The circumcircle centre is returned in (xc,yc) and the radius r
   NOTE: A point on the edge is inside the circumcircle
*/

static bool sf_cram_circum_circle (double xp, double yp, double x1, double y1,
                                   double x2, double y2, double x3, double y3,
                                   double *xc, double *yc, double *rsqr) {
    double m1, m2, mx1, mx2, my1, my2;
    double dx, dy, drsqr;
    double fabsy1y2 = fabs (y1-y2);
    double fabsy2y3 = fabs (y2-y3);

    /* Check for coincident points */
    if (fabsy1y2 < TRI_EPS && fabsy2y3 < TRI_EPS)
        return true;

    if (fabsy1y2 < TRI_EPS) {
        m2 = -(x3-x2)/(y3-y2);
        mx2 = (x2 + x3)/2.0;
        my2 = (y2 + y3)/2.0;
        *xc = (x2 + x1)/2.0;
        *yc = m2*(*xc - mx2) + my2;
    } else if (fabsy2y3 < TRI_EPS) {
        m1 = -(x2-x1)/(y2-y1);
        mx1 = (x1 + x2)/2.0;
        my1 = (y1 + y2)/2.0;
        *xc = (x3 + x2)/2.0;
        *yc = m1*(*xc - mx1) + my1;
    } else {
       m1 = -(x2-x1)/(y2-y1);
       m2 = -(x3-x2)/(y3-y2);
       mx1 = (x1 + x2)/2.0;
       mx2 = (x2 + x3)/2.0;
       my1 = (y1 + y2)/2.0;
       my2 = (y2 + y3)/2.0;
       *xc = (m1*mx1 - m2*mx2 + my2 - my1)/(m1 - m2);
       if (fabsy1y2 > fabsy2y3) {
           *yc = m1*(*xc - mx1) + my1;
       } else {
           *yc = m2*(*xc - mx2) + my2;
       }
    }

    dx = x2 - *xc;
    dy = y2 - *yc;
    *rsqr = dx*dx + dy*dy;

    dx = xp - *xc;
    dy = yp - *yc;
    drsqr = dx*dx + dy*dy;

    return (drsqr - *rsqr) <= TRI_EPS ? true : false;
}

/*
    Takes array xy[np+3] as input (x_i = xy[i*st], y_i = xy[i*st+1],
    returns a list of ntr vertices[3*np] for all detected triangles.
    The vertex array xy must be sorted in the order of increasing x values.
*/

int sf_cram_triangulate (int np, int st, float *xy, int *vertices, int *ntr)
/*< 2-D Delaunay triangulation of an raay of points >*/
{
    bool *complete = NULL;
    int *edges = NULL;
    int nedge = 0;
    int trimax, emax = 200;

    bool inside;
    int i, j, k;
    double xp, yp, x1, y1, x2, y2, x3, y3, xc = 0, yc = 0, r = 0;
    double xmin, xmax, ymin, ymax, xmid, ymid;
    double dx, dy, dmax;

    /* Allocate memory for the completeness list, flag for each triangle */
    trimax = 4*np;
    if ((complete = (bool*)malloc (trimax*sizeof(bool))) == NULL) {
        return 1;
    }

    /* Allocate memory for the edge list */
    if ((edges = (int*)malloc (emax*2*sizeof(int))) == NULL) {
        free (complete);
        return 2;
    }

    /*
       Find the maximum and minimum vertex bounds.
       This is to allow calculation of the bounding triangle
    */
    xmin = xy[0];
    ymin = xy[1];
    xmax = xmin;
    ymax = ymin;
    for (i = 1; i < np; i++) {
        if (xy[i*st] < xmin) xmin = xy[i*st];
        if (xy[i*st] > xmax) xmax = xy[i*st];
        if (xy[i*st + 1] < ymin) ymin = xy[i*st + 1];
        if (xy[i*st + 1] > ymax) ymax = xy[i*st + 1];
    }
    dx = xmax - xmin;
    dy = ymax - ymin;
    dmax = (dx > dy) ? dx : dy;
    xmid = (xmax + xmin)/2.0;
    ymid = (ymax + ymin)/2.0;
    /*
       Set up the supertriangle
       This is a triangle which encompasses all the sample points.
       The supertriangle coordinates are added to the end of the
       vertex list. The supertriangle is the first triangle in
       the triangle list.
    */
    xy[np*st] = xmid - 20*dmax;
    xy[np*st + 1] = ymid - dmax;
    xy[(np + 1)*st] = xmid;
    xy[(np + 1)*st + 1] = ymid + 20*dmax;
    xy[(np + 2)*st] = xmid + 20*dmax;
    xy[(np + 2)*st + 1] = ymid - dmax;
    vertices[0] = np;
    vertices[1] = np + 1;
    vertices[2] = np + 2;
    complete[0] = false;
    *ntr = 1;
    /*
       Include each point one at a time into the existing mesh
    */
    for (i = 0; i < np; i++) {
        xp = xy[i*st];
        yp = xy[i*st + 1];
        nedge = 0;
        /*
           Set up the edge buffer.
           If the point (xp,yp) lies inside the circumcircle then the
           three edges of that triangle are added to the edge buffer
           and that triangle is removed.
        */
        for (j = 0; j < (*ntr); j++) {
            if (complete[j])
                continue;
            x1 = xy[vertices[j*3]*st];
            y1 = xy[vertices[j*3]*st + 1];
            x2 = xy[vertices[j*3 + 1]*st];
            y2 = xy[vertices[j*3 + 1]*st + 1];
            x3 = xy[vertices[j*3 + 2]*st];
            y3 = xy[vertices[j*3 + 2]*st + 1];
            inside = sf_cram_circum_circle (xp, yp, x1, y1, x2, y2, x3, y3, &xc, &yc, &r);
            if (xc < xp && ((xp - xc)*(xp - xc)) > r)
                complete[j] = true;
            if (inside) {
                /* Check that we haven't exceeded the edge list size */
                if (nedge + 3 >= emax) {
                    emax += 100;
                    if ((edges = (int*)realloc (edges, emax*2*sizeof(int))) == NULL) {
                        free (complete);
                        return 3;
                    }
                }
                /* Triangle sides */
                edges[nedge*2] = vertices[j*3];
                edges[nedge*2 + 1] = vertices[j*3 + 1];
                edges[(nedge + 1)*2] = vertices[j*3 + 1];
                edges[(nedge + 1)*2 + 1] = vertices[j*3 + 2];
                edges[(nedge + 2)*2] = vertices[j*3 + 2];
                edges[(nedge + 2)*2 + 1] = vertices[j*3];
                nedge += 3;
                vertices[j*3] = vertices[((*ntr) - 1)*3];
                vertices[j*3 + 1] = vertices[((*ntr) - 1)*3 + 1];
                vertices[j*3 + 2] = vertices[((*ntr) - 1)*3 + 2];
                complete[j] = complete[(*ntr) - 1];
                (*ntr)--;
                j--;
            }
        }
        /*
           Tag multiple edges
           Note: if all triangles are specified anticlockwise then all
                 interior edges are opposite pointing in direction.
        */
        for (j = 0; j < (nedge - 1); j++) {
            for (k = j + 1; k < nedge; k++) {
                if ((edges[j*2] == edges[k*2 + 1]) &&
                    (edges[j*2 + 1] == edges[k*2])) {
                    edges[j*2] = -1;
                    edges[j*2 + 1] = -1;
                    edges[k*2] = -1;
                    edges[k*2 + 1] = -1;
                }
                /* Shouldn't need the following, see note above */
                if ((edges[j*2] == edges[k*2]) &&
                    (edges[j*2 + 1] == edges[k*2 + 1])) {
                    edges[j*2] = -1;
                    edges[j*2 + 1] = -1;
                    edges[k*2] = -1;
                    edges[k*2 + 1] = -1;
                }
            }
        }
        /*
           Form new triangles for the current point
           Skipping over any tagged edges.
           All edges are arranged in clockwise order.
        */
        for (j = 0; j < nedge; j++) {
             if (edges[j*2] < 0 || edges[j*2 + 1] < 0)
                 continue;
             if ((*ntr) >= trimax) {
                 free (complete);
                 free (edges);
                 return 4;;
             }
             vertices[(*ntr)*3] = edges[j*2];
             vertices[(*ntr)*3 + 1] = edges[j*2 + 1];
             vertices[(*ntr)*3 + 2] = i;
             complete[*ntr] = false;
             (*ntr)++;
        }
    }
    /*
       Remove triangles with supertriangle vertices
       These are triangles which have a vertex number greater than np
    */
    for (i = 0; i < (*ntr); i++) {
        if (vertices[i*3] >= np ||
            vertices[i*3 + 1] >= np ||
            vertices[i*3 + 2] >= np) {
            vertices[i*3] = vertices[((*ntr) - 1)*3];
            vertices[i*3 + 1] = vertices[((*ntr) - 1)*3 + 1];
            vertices[i*3 + 2] = vertices[((*ntr) - 1)*3 + 2];
            (*ntr)--;
            i--;
        }
    }

    free (edges);
    free (complete);
    return 0;
}

