/* Data interpolation in 2-D by Thin Plate Spline method. */
/*
  Copyright (C) 2008 University of Texas at Austin

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
/* This TPS implementation is partially based on
   the C++ source code by Jarno Elonen */

#include <float.h>
#include <math.h>

#include <rsf.h>

static float** B;

static void sf_tps_matmult_init (float** bb) {
    B = bb;
}

static void sf_tps_matmult_lop (bool adj, bool add, 
		                int nx, int ny, float* x, float*y) {
    int ix, iy;
    sf_adjnull (adj, add, nx, ny, x, y);
    for (ix = 0; ix < nx; ix++) {
	for (iy = 0; iy < ny; iy++) {
	    if (adj) x[ix] += B[iy][ix] * y[iy];
	    else     y[iy] += B[iy][ix] * x[ix];
	}
    }
}

static double sf_tps_base_func (double r) {
    if (r == 0.0)
        return 0.0;
    else
        return r * r * log (r);
}

static void sf_calc_thin_plate_spline (float **xy, float *d, int nd, float eps, int niter,
                                       float **L, float *x, float *v, float *res, float *grid,
                                       float x0, int nx, float dx, float y0, int ny, float dy) {
    int i, j, k;
    double a = 0.0, elen, h;
    float *xy_i, *xy_j, x_coord, y_coord;

    if (nd < 3)
        return;

/* Fill K (p x p, upper left of L) and calculate
   mean edge length from control points.
   K is symmetrical, thus only half of the 
   coefficients has to be calculated. */
    for (i = 0; i < nd; i++) {
        for (j = i + 1; j < nd; j++) {
            xy_i = xy[i];
            xy_j = xy[j];
            elen = sqrt ((xy_i[0] - xy_j[0]) * (xy_i[0] - xy_j[0]) + 
                         (xy_i[1] - xy_j[1]) * (xy_i[1] - xy_j[1]));
            L[i][j] = sf_tps_base_func (elen);
            L[j][i] = L[i][j];
            a += elen * 2;
        }
    }
    a /= (double)(nd * nd);

/* Fill the rest of L */
    for (i = 0; i < nd; i++) {
/* diagonal: reqularization parameters (lambda * a^2) */
        L[i][i] = eps * (a * a);

/* P (p x 3, upper right) */
        L[i][nd] = 1.0;
        L[i][nd + 1] = xy[i][0];
        L[i][nd + 2] = xy[i][1];

/* P transposed (3 x p, bottom left) */
        L[nd][i] = 1.0;
        L[nd + 1][i] = xy[i][0];
        L[nd + 2][i] = xy[i][1];
    }
/* O (3 x 3, lower right) */
    for (i = nd; i < nd + 3; i++) {
        for (j = nd; j < nd + 3; j++) {
            L[i][j] = 0.0;
        }
    }

/* Fill the right hand vector V */
    for (i = 0; i < nd; i++) {
        v[i] = d[i];
    }
    v[nd] = 0;   
    v[nd + 1] = 0;   
    v[nd + 2] = 0;   

/* Solve the linear system */
    sf_tps_matmult_init (L);
    for (i = 0; i < niter; i++) {
        sf_solver (sf_tps_matmult_lop, sf_cgstep, nd + 3, nd + 3, x, v, i, 
                   "res", res, "end");
        sf_cgstep_close ();
    }

/* Do interpolation */
    for (k = 0; k < ny; k++) {
        for (j = 0; j < nx; j++) {
            x_coord = x0 + j * dx;
            y_coord = y0 + k * dy;
            h = x[nd] + x[nd + 1] * x_coord + x[nd + 2] * y_coord;
            for (i = 0; i < nd; i++) {
                xy_i = xy[i];
                h += x[i] * sf_tps_base_func (sqrt ((xy_i[0] - x_coord) * (xy_i[0] - x_coord) + 
                                                    (xy_i[1] - y_coord) * (xy_i[1] - y_coord)));
            }
            grid[k * nx + j] = h;
        }
    }
}


int main (int argc, char* argv[]) {
    bool verb;
    int id, nk, nd, nm, nt, it, nx, ny, n2, xkey, ykey;
    int niter;
    float *dd, **xy, *hdr, **L, *x, *v, *res, *grid;
    float x0, y0, dx, dy, xmin, xmax, ymin, ymax, f, dt, t0, eps;
    char *header;
    sf_file in, out, head;

    sf_init (argc,argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (!sf_histint (in, "n1", &nd)) sf_error ("Need n1= in in");
    if (nd < 3) sf_error ("Need n1 >= 3");
    if (!sf_histint (in, "n2", &nt)) sf_error ("Need n2= in in");
    if (SF_FLOAT != sf_gettype (in)) sf_error ("Need float input");

    if (!sf_getint("xkey",&xkey)) sf_error("Need xkey=");
    /* x key number */
    if (!sf_getint("ykey",&ykey)) sf_error("Need ykey=");
    /* y key number */

    /* create coordinates */
    xy = sf_floatalloc2 (2, nd);

    header = sf_getstring ("head");
    if (NULL == header) { 
        header = sf_histstring (in, "head");
        if (NULL == header) sf_error ("Need head=");
    }

    head = sf_input (header);

    if (SF_FLOAT != sf_gettype (head)) sf_error ("Need float header");
    if (!sf_histint (head, "n1", &nk)) sf_error ("No n1= in head");
    if (!sf_histint (head, "n2", &n2) || n2 != nd)
        sf_error ("Wrong n2= in head");

    hdr = sf_floatalloc (nk);

    ymin = xmin = +FLT_MAX;
    ymax = xmax = -FLT_MAX;
    /* Determine bounding box and locations of control points */
    for (id = 0; id < nd; id++) {        
        sf_floatread (hdr, nk, head);
        f = hdr[xkey]; 
        if (f < xmin) xmin = f;
        if (f > xmax) xmax = f;
        xy[id][0] = f;
        f = hdr[ykey]; 
        if (f < ymin) ymin = f;
        if (f > ymax) ymax = f;
        xy[id][1] = f;
    }

    sf_fileclose (head);

    /* create model */
    if (!sf_getint ("nx", &nx)) sf_error ("Need nx=");
    /* Number of bins in x */
    if (!sf_getint ("ny", &ny)) sf_error ("Need ny=");
    /* Number of bins in y */

    sf_putint (out, "n1", nx);
    sf_putint (out, "n2", ny);
    sf_putint (out, "n3", nt);

    if (sf_histfloat (in, "o2", &t0)) sf_putfloat (out, "o3", t0);
    if (sf_histfloat (in, "d2", &dt)) sf_putfloat (out, "d3", dt);

    /* let user overwrite */
    sf_getfloat ("xmin", &xmin);
    sf_getfloat ("xmax", &xmax);
    sf_getfloat ("ymin", &ymin);
    /* Grid dimensions */
    sf_getfloat ("ymax", &ymax);

    if (xmax <= xmin) sf_error ("xmax=%f <= xmin=%f", xmax, xmin);
    if (ymax <= ymin) sf_error ("ymax=%f <= ymin=%f", xmax, xmin);

    if (!sf_getfloat ("x0", &x0)) x0 = xmin; 
    if (!sf_getfloat ("y0", &y0)) y0 = ymin; 
    /* grid origin */

    sf_putfloat (out, "o1", x0);
    sf_putfloat (out, "o2", y0);

    if (!sf_getfloat ("dx", &dx)) {
        /* bin size in x */
        if (1 >= nx) sf_error ("Need dx=");
        dx = (xmax - xmin) / (nx - 1);
    }

    if (!sf_getfloat ("dy", &dy)) {
        /* bin size in y */
        if (1 >= nx) {
            dy = dx;
        } else {
            dy = (ymax - ymin) / (ny - 1);
        }
    }

    sf_putfloat (out, "d1", dx);
    sf_putfloat (out, "d2", dy);

    nm = nx * ny;
    grid = sf_floatalloc (nm);
    dd = sf_floatalloc (nd);

    if (!sf_getint ("niter", &niter)) niter = 100;
    /* number of iterations */

    if (!sf_getfloat ("eps", &eps)) eps = 1./ (float)nd;
    /* regularization parameter */
    if (!sf_getbool("verb",&verb)) verb = true;
    /* verbosity flag */

    /* space for interpolation operator and opther part of the system */
    L = sf_floatalloc2 (nd + 3, nd + 3);
    x = sf_floatalloc (nd + 3);
    v = sf_floatalloc (nd + 3);
    res = sf_floatalloc (nd + 3);

    for (it = 0; it < nt; it++) { /* loop over time slices */
        sf_warning ("slice %d of %d", it + 1, nt);
        sf_floatread (dd, nd, in);

        sf_calc_thin_plate_spline (xy, dd, nd, eps, niter,
                                   L, x, v, res, grid,
                                   x0, nx, dx, y0, ny, dy);
        sf_floatwrite (grid, nm, out);
    }


    exit(0);
}
