/* B-spline coefficients for a 3-D (an)isotropic velocity model. */
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

int main (int argc, char* argv[]) {
    int nz, nx, ny, ic, nc;
    float dz, oz, dx, ox, dy, oy;
    size_t i, n, sz;
    float *buf = NULL, *buf2 = NULL;
#ifdef HAVE_SSE
    unsigned char pad[64];
#endif
    sf_file velz, velx = NULL, theta = NULL, phi = NULL, eta = NULL, out;
    bool verb;
    Ugrid z_grid, x_grid, y_grid;
    BCtype_s zBC, xBC, yBC;
    multi_UBspline_3d_s *velspline = NULL;

    sf_init (argc, argv);

    velz = sf_input ("in");
    /* Vertical velocity */
    out = sf_output ("out");
    /* Spline coefficients */

    /* Spatial dimensions */
    if (!sf_histint (velz, "n1", &nz)) sf_error ("No n1= in input");
    if (!sf_histint (velz, "n2", &nx)) sf_error ("No n2= in input");
    if (!sf_histint (velz, "n3", &ny)) sf_error ("No n3= in input");
    if (!sf_histfloat (velz, "d1", &dz)) sf_error ("No d1= in input");
    if (!sf_histfloat (velz, "o1", &oz)) oz = 0.0;
    if (!sf_histfloat (velz, "d2", &dx)) sf_error ("No d2= in input");
    if (!sf_histfloat (velz, "o2", &ox)) ox = 0.0;
    if (!sf_histfloat (velz, "d3", &dy)) sf_error ("No d3= in input");
    if (!sf_histfloat (velz, "o3", &oy)) oy = 0.0;

    if (!sf_getbool ("verb", &verb)) verb = false;
    /* verbosity flag */

    n = (size_t)nz*(size_t)nx*(size_t)ny;
    buf = sf_floatalloc (n);
    nc = 1;

    if (sf_getstring ("vx")) {
        /* Horizontal velocity */
        velx = sf_input ("vx");
        nc++;
    }
    if (sf_getstring ("eta")) {
        /* Anellipticity */
        if (NULL == velx)
            sf_error ("Need vx=, if eta= is given");
        eta = sf_input ("eta");
        nc++;
    } else if (velx)
        sf_error ("Need eta=, if vx= is given");
    if (sf_getstring ("theta")) {
        /* Tilt angle elevation */
        if (NULL == velx)
            sf_error ("Need vx=, if theta= is given");
        theta = sf_input ("theta");
        nc++;
    }
    if (sf_getstring ("phi")) {
        /* Tilt angle azimuth */
        if (NULL == theta)
            sf_error ("Need theta=, if phi= is given");
        phi = sf_input ("phi");
        nc++;
    }

    z_grid.start = oz; z_grid.end = oz + (nz - 1)*dz; z_grid.num = nz;
    x_grid.start = ox; x_grid.end = ox + (nx - 1)*dx; x_grid.num = nx;
    y_grid.start = oy; y_grid.end = oy + (ny - 1)*dy; y_grid.num = ny;
    zBC.lCode = zBC.rCode = NATURAL;
    xBC.lCode = xBC.rCode = NATURAL;
    yBC.lCode = yBC.rCode = NATURAL;
    velspline = create_multi_UBspline_3d_s (y_grid, x_grid, z_grid, yBC, xBC, zBC, nc); 

    /* Read data and compute spline coefficients */
    if (verb)
        sf_warning ("Processing V_z");
    sf_floatread (buf, n, velz);
    if (1 == nc) {
        /* Isotropic case - convert velocity to slowness */
        if (verb)
            sf_warning ("Converting to slowness for isotropic case");
        for (i = 0; i < n; i++)
            buf[i] = 1.0/buf[i];
    } else {
        /* Convert to V_z^2 */
        for (i = 0; i < n; i++)
            buf[i] *= buf[i];
    }
    ic = 0;
    set_multi_UBspline_3d_s (velspline, ic, buf);
    ic++;
    if (velx) {
        if (verb)
            sf_warning ("Processing V_x");
        buf2 = sf_floatalloc (n);
        sf_floatread (buf2, n, velx);
        sf_fileclose (velx);
        /* Convert to V_x^2 */
        for (i = 0; i < n; i++)
            buf2[i] *= buf2[i];
        set_multi_UBspline_3d_s (velspline, ic, buf2);
        ic++;
        /* Convert to (V_z*V_x)^2 */
        for (i = 0; i < n; i++)
            buf[i] *= buf2[i];
    }
    if (eta) {
        if (verb)
            sf_warning ("Processing Eta");
        sf_floatread (buf2, n, eta);
        sf_fileclose (eta);
        /* Convert to -8*eta/(1 + 2*eta)*(V_z*V_x)^2 */
        for (i = 0; i < n; i++) {
            buf2[i] = -8.0*buf2[i]/(1.0 + 2.0*buf2[i]);
            buf2[i] *= buf[i];
        }
        set_multi_UBspline_3d_s (velspline, ic, buf2);
        ic++;
    }
    if (theta) {
        if (verb)
            sf_warning ("Processing Theta");
        sf_floatread (buf, n, theta);
        sf_fileclose (theta);
        /* Convert to radians */
        for (i = 0; i < n; i++)
            buf[i] = buf[i]*SF_PI/180.0;
        set_multi_UBspline_3d_s (velspline, ic, buf);
        ic++;
    }
    if (phi) {
        if (verb)
            sf_warning ("Processing Phi");
        sf_floatread (buf, n, phi);
        sf_fileclose (phi);
        /* Convert to radians */
        for (i = 0; i < n; i++)
            buf[i] = buf[i]*SF_PI/180.0;
        set_multi_UBspline_3d_s (velspline, ic, buf);
        ic++;
    }

    if (buf2)
        free (buf2);
    free (buf);

    sz = (size_t)sizeof(multi_UBspline_3d_s) +
         (size_t)velspline->nc;
#ifdef HAVE_SSE
    if (sizeof(multi_UBspline_3d_s) % 64)
        sz += 64 - (sizeof(multi_UBspline_3d_s) % 64);
#endif
    /* Make output a 1-D file of coefficients */
    sf_unshiftdim2 (velz, out, 1);
    /* Set up output */
    sf_settype (out, SF_UCHAR);
    sf_putlargeint (out, "n1", sz);
    sf_putfloat (out, "o1", 0.0);
    sf_putfloat (out, "d1", 1.0);
    sf_putstring (out, "label1", "Spline coefficients");
    sf_putstring (out, "unit1", "");
    sf_putstring (out, "label2", "");
    sf_putstring (out, "unit2", "");

    sf_putint (out, "Nz", nz);
    sf_putfloat (out, "Oz", oz);
    sf_putfloat (out, "Dz", dz);
    sf_putint (out, "Nx", nx);
    sf_putfloat (out, "Ox", ox);
    sf_putfloat (out, "Dx", dx);
    sf_putint (out, "Ny", ny);
    sf_putfloat (out, "Oy", oy);
    sf_putfloat (out, "Dy", dy);

    sf_putint (out, "Nc", nc);
    sf_putstring (out, "splines", "y");

    if (verb) {
        sf_warning ("Number of spline coefficients: %lu",
                    velspline->nc/(size_t)sizeof(float));
        sf_warning ("Writing spline coefficients");
    }

    sf_ucharwrite ((unsigned char*)velspline, (size_t)sizeof(multi_UBspline_3d_s), out);
#ifdef HAVE_SSE
    if (sizeof(multi_UBspline_3d_s) % 64)
        sf_ucharwrite (pad, (size_t)(64 - (sizeof(multi_UBspline_3d_s) % 64)), out);
#endif
    sf_ucharwrite ((unsigned char*)velspline->coefs, (size_t)velspline->nc, out);
    destroy_Bspline (velspline);

    return 0;
}

