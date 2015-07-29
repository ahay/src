/* 2-D angle-domain Kirchhoff migration based on escape tables. */
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

#include "cram_gather2.h"
#include "cram_point2.h"
#include "esc_point2.h"

int main (int argc, char* argv[]) {
    int iz, ix, nz, nx, na, nt, ts, th;
    float dt, da, t0, a0;
    float dz, z0, dx, x0, zd, z, x, zf, vconst = 1.5, smax, hmax;
    float oazmin = 180.0, oazmax = 180.0, dazmin = 180.0, dazmax = 180.0;
    float oaz, daz, ds, dh;
    float **esc;
    sf_file esct, data = NULL, vz = NULL;
    sf_file oimag = NULL, dimag = NULL, fimag = NULL,
            osmap = NULL, dsmap = NULL, oimap = NULL, dimap = NULL;

    bool mute, sqsmb;
    sf_cram_data2 cram_data;
    sf_cram_survey2 cram_survey;
    sf_cram_slowness2 cram_slowness;
    sf_cram_rbranch2 cram_rbranch;
    sf_cram_point2 cram_point;
    sf_cram_gather2 cram_gather;

    sf_init (argc, argv);

    esct = sf_input ("in");
    /* Escape tables */

    /* Phase space dimensions (also, imaging dimensions) */
    if (!sf_histint (esct, "n1", &na)) sf_error ("No n1= in input");
    if (na != ESC2_NUM) sf_error ("Need n1=%d in input", ESC2_NUM);
    if (!sf_histint (esct, "n2", &na)) sf_error ("No n2= in input");
    if (!sf_histint (esct, "n3", &nz)) sf_error ("No n3= in input");
    if (!sf_histint (esct, "n4", &nx)) sf_error ("No n4= in input");
    if (!sf_histfloat (esct, "d2", &da)) sf_error ("No d2= in input");
    if (!sf_histfloat (esct, "o2", &a0)) sf_error ("No o2= in input");
    if (!sf_histfloat (esct, "d3", &dz)) sf_error ("No d3= in input");
    if (!sf_histfloat (esct, "o3", &z0)) sf_error ("No o3= in input");
    if (!sf_histfloat (esct, "d4", &dx)) sf_error ("No d4= in input");
    if (!sf_histfloat (esct, "o4", &x0)) sf_error ("No o4= in input");

    /* Surface depth */
    zd = z0 + 0.25*dz;

    if (!sf_getbool ("mute", &mute)) mute = false;
    /* y - mute signal in constant z plane before stacking */
    if (!sf_getbool ("sqsmb", &sqsmb)) sqsmb = false;
    /* y - output energy traces instead of semblance */

    if (mute) {
        if (!sf_getfloat ("oazmin", &oazmin)) oazmin = 180.0;
        /* Maximum allowed scattering angle at z min */
        if (!sf_getfloat ("oazmax", &oazmax)) oazmax = 180.0;
        /* Maximum allowed scattering angle at z max */
        if (!sf_getfloat ("dazmin", &dazmin)) dazmin = 180.0;
        /* Maximum allowed dip angle (abs.value) at z min */
        if (!sf_getfloat ("dazmax", &dazmax)) dazmax = 180.0;
        /* Maximum allowed dip angle (abs.value) at z max */
        if (oazmin < 0.0) oazmin = 0.0;
        if (oazmax < 0.0) oazmax = 0.0;
        if (dazmin < 0.0) dazmin = 0.0;
        if (dazmax < 0.0) dazmax = 0.0;
    }

    if (!sf_getint ("ts", &ts)) ts = 3;
    /* Tapering length at the edges of the source direction */
    if (!sf_getint ("th", &th)) th = 5;
    /* Tapering length at the edges of the receiver direction */
    
    esc = sf_floatalloc2 (ESC2_NUM, na);

    if (sf_getstring ("data")) {
        /* Processed prestack data */
        data = sf_input ("data");
    } else {
        sf_error ("Need data=");
    }

    /* Data dimensions */
    if (!sf_histint (data, "n1", &nt)) sf_error ("No n1= in data");
    if (!sf_histfloat (data, "d1", &dt)) sf_error ("No d1= in data");
    if (!sf_histfloat (data, "o1", &t0)) sf_error ("No o1= in data");

    if (sf_getstring ("vz")) {
        /* Velocity model for amplitude weights */
        vz = sf_input ("vz");
    } else {
        if (!sf_getfloat ("vconst", &vconst)) vconst = 1.5;
        /* Constant velocity, if vz= is not used */
    }

    /* Data object */
    cram_data = sf_cram_data2_init (data, NULL);
    /* Survey object */
    cram_survey = sf_cram_survey2_init (data);

    ds = sf_cram_survey2_get_src_sampling (cram_survey);
    dh = sf_cram_survey2_get_rcv_sampling (cram_survey);
    if (!sf_getfloat ("smax", &smax)) smax = 10*fabsf (ds);
    /* Maximum allowed width of the shot ray branch  */
    if (!sf_getfloat ("hmax", &hmax)) hmax = 20*fabsf (dh);
    /* Maximum allowed width of the receiver ray branch  */

    /* Slowness object */
    cram_slowness = sf_cram_slowness2_init (vz, vconst);
    /* Exit ray branches object */
    cram_rbranch = sf_cram_rbranch2_init (na, zd, t0 + (nt - 1)*dt, cram_slowness);
    sf_cram_rbranch2_set_maxwidth (cram_rbranch, smax > hmax ? hmax : smax);
    /* Subsurface point image object */
    cram_point = sf_cram_point2_init (na, a0, da, cram_data, cram_survey, cram_slowness,
                                      cram_rbranch);
    sf_cram_point2_set_taper (cram_point, ts, th);
    sf_cram_point2_set_shmax (cram_point, smax, hmax);
    /* Image and gathers accumulator object */
    cram_gather = sf_cram_gather2_init (na, nz, a0, da, oazmax, dazmax);
    sf_cram_gather2_set_sqsmb (cram_gather, sqsmb);

    oimag = sf_output ("out");
    /* Scattering angle gathers (angle, z, x) */
    sf_cram_gather2_setup_oangle_output (cram_gather, esct, oimag, false);
    if (sf_getstring ("imap")) {
        /* Scattering gathers illumination (angle, z, x) */
        oimap = sf_output ("imap");
        sf_cram_gather2_setup_oangle_output (cram_gather, esct, oimap, true);
    }
    if (sf_getstring ("smap")) {
        /* Scattering gathers semblance (angle, z, x) */
        osmap = sf_output ("smap");
        sf_cram_gather2_setup_oangle_output (cram_gather, esct, osmap, false);
    }

    if (sf_getstring ("dipagath")) {
        /* Dip angle gathers (angle, z, x) */
        dimag = sf_output ("dipagath");
        sf_cram_gather2_setup_dangle_output (cram_gather, esct, dimag, false);
        if (sf_getstring ("dipimap")) {
            /* Dip gathers illumination (angle, z, x) */
            dimap = sf_output ("dipimap");
            sf_cram_gather2_setup_dangle_output (cram_gather, esct, dimap, true);
        }
        if (sf_getstring ("dipsmap")) {
            /* Dip gathers semblance (angle, z, x) */
            dsmap = sf_output ("dipsmap");
            sf_cram_gather2_setup_dangle_output (cram_gather, esct, dsmap, false);
        }
    }

    if (sf_getstring ("full")) {
        /* Full image (scattering angle, dip angle, z, x) */
        fimag = sf_output ("full");
        sf_cram_gather2_setup_full_output (cram_gather, fimag, false);
    }

    for (ix = 0; ix < nx; ix++) { /* Loop over image x */
        x = x0 + ix*dx;
        sf_warning ("Lateral %d of %d (%g)", ix + 1, nx, x);
        for (iz = 0; iz < nz; iz++) { /* Loop over image z */
            z = z0 + iz*dz;
            /* Escape variables for this (z,x) location */
            sf_floatread (esc[0], ESC2_NUM*na, esct);
            if (mute) {
                zf = (nz != 1) ? iz/(float)(nz - 1) : 0.0;
                /* Maximum scattering angle for this z */
                oaz = oazmin + zf*(oazmax - oazmin);
                /* Maximum dip angle for this z */
                daz = dazmin + zf*(dazmax - dazmin);
                sf_cram_point2_set_mute (cram_point, oaz, daz);
            }
            /* Image for one (z,x) location */
            sf_cram_point2_compute (cram_point, z, x, esc);
            /* Add to the angle gathers */
            sf_cram_gather2_set_point (cram_gather, iz, cram_point);
            if (fimag)
                sf_cram_gather2_full_output (cram_gather, cram_point, fimag, NULL);
        } /* iz */
        if (dimag) {
            sf_warning ("Lateral %g: saving dip gathers", x0 + ix*dx);
            sf_cram_gather2_dangle_output (cram_gather, dimag, dsmap, dimap);
        }
        sf_warning ("Lateral %g: saving angle gathers", x0 + ix*dx);
        sf_cram_gather2_oangle_output (cram_gather, oimag, osmap, oimap);
    } /* ix */

    sf_cram_gather2_close (cram_gather);
    sf_cram_point2_close (cram_point);
    sf_cram_rbranch2_close (cram_rbranch);
    sf_cram_slowness2_close (cram_slowness);
    sf_cram_survey2_close (cram_survey);
    sf_cram_data2_close (cram_data);

    sf_fileclose (data);
    if (vz)
        sf_fileclose (vz);
    if (dimag)
        sf_fileclose (dimag);
    if (fimag)
        sf_fileclose (fimag);
    if (osmap)
        sf_fileclose (osmap);
    if (dsmap)
        sf_fileclose (dsmap);
    if (oimap)
        sf_fileclose (oimap);
    if (dimap)
        sf_fileclose (dimap);

    free (esc[0]);
    free (esc);

    return 0;
}

