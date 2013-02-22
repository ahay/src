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

#include "esc_bgrid3.h"
#include "esc_tracer3.h"

static char *dbname = "scblock";

int main (int argc, char* argv[]) {
    size_t i, n, nc, nnc;
    int nz, nx, ny, nsb, nsa, ic, ia, iscz, iscx, iscy;
    int nsz, nsx, nsy, nscz, nscx, nscy;
    float dz, oz, dx, ox, dy, oy;
    float oscz, oscx, oscy, dscz, dscx, dscy;
    float zmin, zmax, xmin, xmax, ymin, ymax;
    multi_UBspline_3d_s **splines[ESC3_SIDE_NUM];
    sf_file scmap, vspline = NULL, scblock = NULL, out;

    char sbuf[ESC3_BSC_BMAX];
    char *bname = NULL;
    bool verb, parab;
    sf_esc_slowness3 esc_slow;
    sf_esc_tracer3 esc_tracer;
    sf_esc_bgrid3 esc_bgrid;
    sf_simtab simtab;

    sf_init (argc, argv);

    scmap = sf_input ("in");
    /* Map of supercells to build */
    out = sf_output ("out");
    /* List of supercell objects */

    if (!sf_getbool ("verb", &verb)) verb = false;
    /* verbosity flag */

    bname = sf_getstring ("bname");
    /* Prefix for supercell block file */
    if (NULL == bname)
        bname = dbname;

    if (!sf_histint (scmap, "n1", &nscz)) sf_error ("No n1= in input");
    /* Number of supercells in z */
    if (!sf_histint (scmap, "n2", &nscx)) sf_error ("No n2= in input");
    /* Number of supercells in x */
    if (!sf_histint (scmap, "n3", &nscy)) sf_error ("No n3= in input");
    /* Number of supercells in y */

    if (!sf_histfloat (scmap, "Dscz", &dscz)) sf_error ("No Dscz= in input");
    /* Supercell size in z  */
    if (!sf_histfloat (scmap, "Dscx", &dscx)) sf_error ("No Dscx= in input");
    /* Supercell size in x */
    if (!sf_histfloat (scmap, "Dscy", &dscy)) sf_error ("No Dscy= in input");
    /* Supercell size in y */
    if (!sf_histfloat (scmap, "o1", &oscz)) sf_error ("No Dscz= in input");
    /* Supercell start in z  */
    if (!sf_histfloat (scmap, "o2", &oscx)) sf_error ("No Dscx= in input");
    /* Supercell start in x */
    if (!sf_histfloat (scmap, "o3", &oscy)) sf_error ("No Dscy= in input");
    /* Supercell start in y */
    if (!sf_histint (scmap, "Nz", &nsz)) sf_error ("No Nz= in input");
    /* Supercell sampling in z */
    if (!sf_histint (scmap, "Nx", &nsx)) sf_error ("No Nx= in input");
    /* Supercell sampling in x */
    if (!sf_histint (scmap, "Ny", &nsy)) sf_error ("No Ny= in input");
    /* Supercell sampling in y */
    if (!sf_histint (scmap, "Na", &nsa)) sf_error ("No Na= in input");
    /* Supercell sampling in a */
    if (!sf_histint (scmap, "Nb", &nsb)) sf_error ("No Nb= in input");
    /* Supercell sampling in b */

    if (!sf_getstring ("vspl")) sf_error ("Need vspl=");
    /* Spline coefficients for velocity model */
    vspline = sf_input ("vspl");

    /* Sampling size in the supercell */
    dz = dscz/(float)nsz;
    dx = dscx/(float)nsx;
    dy = dscy/(float)nsy;

    /* Slowness components module [(an)isotropic] */
    esc_slow = sf_esc_slowness3_init (vspline, verb);

    /* Global limits */
    zmin = sf_esc_slowness3_oz (esc_slow);
    zmax = sf_esc_slowness3_oz (esc_slow) +
           (sf_esc_slowness3_nz (esc_slow) - 1)*
           sf_esc_slowness3_dz (esc_slow);
    xmin = sf_esc_slowness3_ox (esc_slow);
    xmax = sf_esc_slowness3_ox (esc_slow) +
           (sf_esc_slowness3_nx (esc_slow) - 1)*
           sf_esc_slowness3_dx (esc_slow);
    ymin = sf_esc_slowness3_oy (esc_slow);
    ymax = sf_esc_slowness3_oy (esc_slow) +
           (sf_esc_slowness3_ny (esc_slow) - 1)*
           sf_esc_slowness3_dy (esc_slow);

    if (!sf_getbool ("parab", &parab)) parab = true;
    /* y - use parabolic approximation of trajectories, n - straight line */

    /* Ray tracer object */
    esc_tracer = sf_esc_tracer3_init (esc_slow, NULL, 0.0, NULL);
    sf_esc_tracer3_set_parab (esc_tracer, parab);

    sf_shiftdim (scmap, out, 1);
    sf_settype (out, SF_UCHAR);
    sf_putint (out, "n1", ESC3_BSC_BMAX);
    sf_putfloat (out, "o1", 0.0);
    sf_putfloat (out, "d1", 1.0);
    sf_putstring (out, "label1", "Array of pointers to supercells");
    sf_putstring (out, "unit1", "");

    simtab = sf_getpars ();

    n = (size_t)nscy*(size_t)nscx*(size_t)nscz;
    i = 0;
    nnc = 0;
    /* Loop over supercells */
    for (iscy = 0; iscy < nscy; iscy++) {
        oy = ymin + oscy + iscy*dscy;
        ny = nsy + 1;
        while ((oy + (ny - 1.01)*dy) > ymax)
            ny--;
        for (iscx = 0; iscx < nscx; iscx++) {
            ox = xmin + oscx + iscx*dscx;
            nx = nsx + 1;
            while ((ox + (nx - 1.01)*dx) > xmax)
                nx--;
            for (iscz = 0; iscz < nscz; iscz++) {
                oz = zmin + oscz + iscz*dscz;
                nz = nsz + 1;
                while ((oz + (nz - 1.01)*dz) > zmax)
                    nz--;
                if (verb)
                    sf_warning ("Processing supercell block %lu of %lu at y=%g, x=%g, z=%g [%dx%dx%dx%dx%d]",
                                i + 1, n, oy, ox, oz, nz, nx, ny, nsb, nsa);
                esc_bgrid = sf_esc_bgrid3_init (nz, nx, ny, nsa, nsb,
                                                oz, ox, oy, dz, dx, dy,
                                                esc_tracer);
                sf_esc_bgrid3_set_verb (esc_bgrid, verb);
                /* Compute escape function on all of the six faces of the
                   supercell */
                splines[ESC3_SIDE_TOP] = sf_esc_bgrid3_compute_topbottom (esc_bgrid, true);
                splines[ESC3_SIDE_BOTTOM] = sf_esc_bgrid3_compute_topbottom (esc_bgrid, false);
                splines[ESC3_SIDE_LEFT] = sf_esc_bgrid3_compute_leftright (esc_bgrid, true);
                splines[ESC3_SIDE_RIGHT] = sf_esc_bgrid3_compute_leftright (esc_bgrid, false);
                splines[ESC3_SIDE_NEAR] = sf_esc_bgrid3_compute_nearfar (esc_bgrid, true);
                splines[ESC3_SIDE_FAR] = sf_esc_bgrid3_compute_nearfar (esc_bgrid, false);
                nc = 0;
                /* Compute total size occupied by the spline function */
                for (ic = 0; ic < ESC3_SIDE_NUM; ic++) {
                    for (ia = 0; ia < nsa; ia++) {
                        nc += (size_t)sizeof(multi_UBspline_3d_s) +
                              (size_t)splines[ic][ia]->nc;
                        nnc += nc;
                    }
                }
                /* Create a separate output file for each supercell */
                snprintf (sbuf, ESC3_BSC_BMAX, "%s_%f_%f_%f.rsf", bname, oz, ox, oy);
                sbuf[ESC3_BSC_BMAX - 1] = '\0';
                if (verb)
                    sf_warning ("Writing supercell block splines, %g Mb", nc*1e-6);
                /* Add a key for the new filename - the only way to
                   create a new RSF file without a command line parameter */
                sf_simtab_enter (simtab, bname, sbuf);
                scblock = sf_output (bname);
                sf_settype (scblock, SF_UCHAR);
                sf_putlargeint (scblock, "n1", nc);
                sf_putfloat (scblock, "o1", 0.0);
                sf_putfloat (scblock, "d1", 1.0);
                sf_putstring (scblock, "label1", "Spline coefficients");
                sf_putstring (scblock, "unit1", "");
                sf_putint (scblock, "n2", 1);
                sf_putint (scblock, "n3", 1);
                sf_putint (scblock, "n4", 1);
                sf_putint (scblock, "Nz", nz);
                sf_putint (scblock, "Nx", nx);
                sf_putint (scblock, "Ny", ny);
                sf_putint (scblock, "Na", nsa);
                sf_putint (scblock, "Nb", nsb);
                sf_putfloat (scblock, "Oz", oz);
                sf_putfloat (scblock, "Ox", ox);
                sf_putfloat (scblock, "Oy", oy);
                sf_putfloat (scblock, "Dz", dz);
                sf_putfloat (scblock, "Dx", dx);
                sf_putfloat (scblock, "Dy", dy);
                sf_putlargeint (scblock, "Nc", nc);
                for (ic = 0; ic < ESC3_SIDE_NUM; ic++) {
                    for (ia = 0; ia < nsa; ia++) {
                        sf_ucharwrite ((unsigned char*)splines[ic][ia],
                                       (size_t)sizeof(multi_UBspline_3d_s), scblock);
                    }
                }
                for (ic = 0; ic < ESC3_SIDE_NUM; ic++) {
                    for (ia = 0; ia < nsa; ia++) {
                        sf_ucharwrite ((unsigned char*)splines[ic][ia]->coefs,
                                       (size_t)splines[ic][ia]->nc, scblock);
                    }
                }
                /* Done with this cell */
                sf_fileclose (scblock);
                sf_ucharwrite ((unsigned char*)sbuf, ESC3_BSC_BMAX, out);
                sf_esc_bgrid3_close (esc_bgrid);
                /* Delete splines */
                for (ic = 0; ic < ESC3_SIDE_NUM; ic++) {
                    for (ia = 0; ia < nsa; ia++) {
                        destroy_Bspline (splines[ic][ia]);
                    }
                    free (splines[ic]);
                }
                i++;
            }
        }
    }
    if (verb)
        sf_warning ("%lu supercell block processed, %g Mb of spline coefficients produced", i, nnc*1e-6);

    sf_esc_tracer3_close (esc_tracer);
    sf_esc_slowness3_close (esc_slow);

    sf_fileclose (vspline);

    return 0;
}

