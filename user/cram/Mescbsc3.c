/* Prepare supercells for stitching escape tables in 3-D. */
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

#include "einspline.h"
#include "esc_helper.h"
#include "esc_point3.h"
#include "esc_tracer3.h"

int main (int argc, char* argv[]) {
    size_t nc = 0;
    int nz, nx, ny, nb, na, nab, iz, ix, iy, ib, ia, iab, i;
    float dz, oz, dx, ox, dy, oy, da, oa, db, ob, oab;
    float z, x, y, a, b, ae, be, vf[3], vt[3], q[4], md;
    float ****e;
    Ugrid z_grid, x_grid, y_grid;
    BCtype_s zBC, xBC, yBC;
    multi_UBspline_3d_s *zxyspline;
    sf_file adom, vspline = NULL, out;

    bool verb, parab;
    sf_esc_slowness3 esc_slow;
    sf_esc_tracer3 esc_tracer;
    sf_esc_point3 esc_point;

    sf_init (argc, argv);

    if (!sf_stdin ()) {
        adom = NULL;
    } else {
        adom = sf_input ("in");
        /* Angle (a*b) domain */
    }

    out = sf_output ("out");
    /* Spline coefficients for local escape functions in supercells */

    if (!sf_getbool ("verb", &verb)) verb = false;
    /* verbosity flag */

    if (!sf_getint ("nz", &nz)) sf_error ("Need nz=");
    /* Number of samples in z axis */
    if (!sf_getfloat ("oz", &oz)) sf_error ("Need oz=");
    /* Beginning of z axis */
    if (!sf_getfloat ("dz", &dz)) sf_error ("Need dz=");
    /* Sampling of z axis */
    if (!sf_getint ("nx", &nx)) sf_error ("Need nx=");
    /* Number of samples in x axis */
    if (!sf_getfloat ("ox", &ox)) sf_error ("Need ox=");
    /* Beginning of x axis */
    if (!sf_getfloat ("dx", &dx)) sf_error ("Need dx=");
    /* Sampling of x axis */
    if (!sf_getint ("ny", &ny)) sf_error ("Need ny=");
    /* Number of samples in y axis */
    if (!sf_getfloat ("oy", &oy)) sf_error ("Need oy=");
    /* Beginning of y axis */
    if (!sf_getfloat ("dy", &dy)) sf_error ("Need dy=");
    /* Sampling of y axis */

    if (adom) {
        if (!sf_histint (adom, "Nb", &nb)) sf_error ("No Na= in input");
        if (!sf_histint (adom, "Na", &na)) sf_error ("No Nb= in input");
        if (!sf_histint (adom, "n1", &nab)) sf_error ("No n1= in input");
        if (!sf_histfloat (adom, "o1", &oab)) sf_error ("No o1= in input");
    } else {
        if (!sf_getint ("na", &na)) sf_error ("Need na=");
        /* Number of samples in azimuth dimension */
        if (!sf_getint ("nb", &nb)) sf_error ("Need nb=");
        /* Number of samples in inclination dimension */
        nab = na*nb;
        oab = 0.0;
    }
    da = 360.0/(float)na;
    /* Sampling of azimtuh dimension */
    oa = 0.5*da;
    /* Beginning of azimuth dimension */
    db = 180.0/(float)nb;
    /* Sampling of inclination dimension */
    ob = 0.5*db;
    /* Beginning of inclination dimension */

    if (!sf_getstring ("vspl")) sf_error ("Need vspl=");
    /* Spline coefficients for velocity model */
    vspline = sf_input ("vspl");

    /* Slowness components module [(an)isotropic] */
    esc_slow = sf_esc_slowness3_init (vspline, verb);

    if (!sf_getbool ("parab", &parab)) parab = true;
    /* y - use parabolic approximation of trajectories, n - straight line */

    /* Ray tracer object */
    esc_tracer = sf_esc_tracer3_init (esc_slow, NULL, 0.0, NULL);
    sf_esc_tracer3_set_parab (esc_tracer, parab);

    /* Set up output */
    sf_settype (out, SF_UCHAR);
    sf_putfloat (out, "o1", 0.0);
    sf_putfloat (out, "d1", 1.0);
    sf_putstring (out, "label1", "Spline coefficients");
    sf_putstring (out, "unit1", "");
    sf_putint (out, "n2", nab);
    sf_putfloat (out, "d2", 1.0);
    sf_putfloat (out, "o2", oab);
    sf_putstring (out, "label2", "Inclination*Azimuth");
    sf_putstring (out, "unit2", "");

    sf_putint (out, "Nz", nz);
    sf_putfloat (out, "Dz", dz);
    sf_putfloat (out, "Oz", oz);
    sf_putint (out, "Nx", nx);
    sf_putfloat (out, "Dx", dx);
    sf_putfloat (out, "Ox", ox);
    sf_putint (out, "Ny", ny);
    sf_putfloat (out, "Dy", dy);
    sf_putfloat (out, "Oy", oy);
    sf_putint (out, "Na", na);
    sf_putfloat (out, "Da", da);
    sf_putfloat (out, "Oa", oa);
    sf_putint (out, "Nb", nb);
    sf_putfloat (out, "Db", db);
    sf_putfloat (out, "Ob", ob);

    if (verb) {
        sf_warning ("Spatial domain dimensions: nz=%d, z=[%g, %g]", nz,
                    oz, oz + (nz - 1)*dz);
        sf_warning ("Spatial domain dimensions: nx=%d, x=[%g, %g]", nx,
                    ox, ox + (nx - 1)*dx);
        sf_warning ("Spatial domain dimensions: ny=%d, y=[%g, %g]", ny,
                    oy, oy + (ny - 1)*dy);
        ia = ((int)(oab + 0.5)) / na;
        ib = ((int)(oab + 0.5)) % na;
        sf_warning ("Angular domain starts at: a=%g, b=%g",
                    oa + ia*da, ob + ib*db);
        ia = ((int)(oab + nab - 0.5)) / na;
        ib = ((int)(oab + nab - 0.5)) % na;
        sf_warning ("Angular domain ends at: a=%g, b=%g (%d angles total)",
                    oa + ia*da, ob + ib*db, nab);
        sf_warning ("Angular domain sampling: da=%g, db=%g", da, db);
    }

    da *= SF_PI/180.0;
    oa *= SF_PI/180.0;
    db *= SF_PI/180.0;
    ob *= SF_PI/180.0;

    z_grid.start = oz;
    z_grid.end = oz + (nz - 1)*dz;
    z_grid.num = nz;
    x_grid.start = ox;
    x_grid.end = ox + (nx - 1)*dx;
    x_grid.num = nx;
    y_grid.start = oy;
    y_grid.end = oy + (ny - 1)*dy;
    y_grid.num = ny;
    zBC.lCode = zBC.rCode = NATURAL;
    xBC.lCode = xBC.rCode = NATURAL;
    yBC.lCode = yBC.rCode = NATURAL;

    if (!sf_getfloat ("md", &md)) md = dz;
    /* Half-width of a supercell */
    sf_esc_tracer3_set_mdist (esc_tracer, md);

    sf_putfloat (out, "Mdist", md);

    e = sf_floatalloc4 (nz, nx, ny, ESC3_NUM + 4);

    esc_point = sf_esc_point3_init ();

    /* Loop over angle domain and create constant-angle
       local escape solutions */
    for (iab = 0; iab < nab; iab++) {
        ia = (int)(oab + iab + 0.5) / na;
        ib = (int)(oab + iab + 0.5) % na;
        a = oa + ia*da;
        b = ob + ib*db;
        /* Initial phase vector */
        vf[0] = cosf (b);
        vf[1] = sinf (b)*cosf (a);
        vf[2] = sinf (b)*sinf (a);
        for (iy = 0; iy < ny; iy++) {
            y = oy + iy*dy;
            sf_esc_tracer3_set_ymin (esc_tracer, y - md);
            sf_esc_tracer3_set_ymax (esc_tracer, y + md);
            for (ix = 0; ix < nx; ix++) {
                x = ox + ix*dx;
                sf_esc_tracer3_set_xmin (esc_tracer, x - md);
                sf_esc_tracer3_set_xmax (esc_tracer, x + md);
                for (iz = 0; iz < nz; iz++) {
                    z = oz + iz*dz;
                    sf_esc_tracer3_set_zmin (esc_tracer, z - md);
                    sf_esc_tracer3_set_zmax (esc_tracer, z + md);
                    sf_esc_tracer3_compute (esc_tracer, z, x, y, b, a,
                                            0.0, 0.0, esc_point, &ae, &be);
                    /* Copy escape values to the output buffer */
                    for (i = 0; i < ESC3_NUM; i++)
                        e[i][iy][ix][iz] = sf_esc_point3_get_esc_var (esc_point, i);
                    e[ESC3_Z][iy][ix][iz] -= z;
                    e[ESC3_X][iy][ix][iz] -= x;
                    e[ESC3_Y][iy][ix][iz] -= y;
                    /* New phase vector */
                    vt[0] = cosf (be);
                    vt[1] = sinf (be)*cosf (ae);
                    vt[2] = sinf (be)*sinf (ae);
                    /* Encode rotation to the new phase vector */
                    sf_quat_vecrot (vf, vt, q);
                    e[ESC3_NUM][iy][ix][iz] = q[0];
                    e[ESC3_NUM + 1][iy][ix][iz] = q[1];
                    e[ESC3_NUM + 2][iy][ix][iz] = q[2];
                    e[ESC3_NUM + 3][iy][ix][iz] = q[3];
                } /* z */
            } /* x */
        } /* y */
        /* Create for this constant-azimuth volume */
        zxyspline = create_multi_UBspline_3d_s (y_grid, x_grid, z_grid, yBC, xBC, zBC, ESC3_NUM + 4);
        for (i = 0; i < (ESC3_NUM + 4); i++) {
            set_multi_UBspline_3d_s (zxyspline, i, &e[i][0][0][0]);
        }
        if (0 == iab) {
            sf_putint (out, "n1", (size_t)sizeof(multi_UBspline_3d_s) +
                                  (size_t)zxyspline->nc);
        }
        sf_ucharwrite ((unsigned char*)zxyspline,
                       (size_t)sizeof(multi_UBspline_3d_s), out);
        sf_ucharwrite ((unsigned char*)zxyspline->coefs,
                       (size_t)zxyspline->nc, out);
        fflush (sf_filestream (out));
        nc += (size_t)zxyspline->nc;
        destroy_Bspline (zxyspline);
    } /* a*b */
    if (verb) {
        sf_warning (".");
        sf_warning ("%lu blocks computed, %g Mb of spline coefficients produced",
                    (size_t)nab, nc*1e-6);
    }

    free (e[0][0][0]);
    free (e[0][0]);
    free (e[0]);
    free (e);

    sf_esc_point3_close (esc_point);
    sf_esc_tracer3_close (esc_tracer);
    sf_esc_slowness3_close (esc_slow);

    sf_fileclose (vspline);

    return 0;
}

