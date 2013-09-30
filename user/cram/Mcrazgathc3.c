/* Collapse/stack (partially) azimuthal axis of 3-D angle-domain migration angle gathers. */
/*
  Copyright (C) 2013 University of Texas at Austin

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

int main (int argc, char* argv[]) {
    int i, j, n, nd, na, nz, nth, wd, ia, iz, iia, iiz, z0, z1;
    float oa, dz;
    float **gatherin, **gatherout;
    sf_file gath, rgath;

    bool verb, dip, norm;

    sf_init (argc, argv);

    gath = sf_input ("in");
    /* Angle gathers */

    if (!sf_histint (gath, "n1", &na)) sf_error ("No n1= in input");
    /* Angle */
    if (!sf_histint (gath, "n2", &nz)) sf_error ("No n2= in input");
    /* Azimuth */
    if (!sf_histfloat (gath, "o1", &oa)) sf_error ("No d1= in input");
    if (!sf_histfloat (gath, "d2", &dz)) sf_error ("No d2= in input");

    dip = oa < 0; /* Dip or scattering angle gathers */

    if (!sf_histint (gath, "n3", &nd)) sf_error ("No n3= in input");
    /* Depth for counting gathers */
    n = sf_leftsize (gath, 2);

    if (!sf_getbool ("verb", &verb)) verb = true;
    /* verbosity flag */
    if (!sf_getbool ("norm", &norm)) norm = false;
    /* y - normalize after stacking */

    if (!sf_getint ("nth", &nth)) nth = 10;
    /* leave every nth azimuth */
    if (!sf_getint ("wd", &wd)) wd = 5;
    /* half-width of stacking base (total base is 2*wd + 1) */

    if (wd < 0)
        wd = 0;
    if (wd > nz/2)
        wd = nz/2;
    if (nth < 1)
        nth = 1;

    rgath = sf_output ("out");
    /* Reduced angle gather */
    sf_putint (rgath, "n2", nz/nth);
    sf_putfloat (rgath, "d2", dz*nth);

    gatherin = sf_floatalloc2 (na, nz);
    gatherout = sf_floatalloc2 (na, nz/nth);

    for (i = 0; i < n; i++) { /* Loop over all depths in all gathers */
        if (verb && 0 == (i % nd))
            sf_warning ("Gather %d of %d", i/nd + 1, n/nd);
        sf_floatread (gatherin[0], na*nz, gath);
        memset (gatherout[0], 0, (nz/nth)*na*sizeof(float));
        for (j = 0; j < nz/nth; j++) { /* Loop over reduced azimuths */
            /* Summation base */
            z0 = j*nth - wd;
            z1 = j*nth + wd;
            for (iz = z0; iz <= z1; iz++) { /* Loop over azimuthal summation base */
                for (ia = 0; ia < na; ia++) { /* Loop over all angles in input */
                    iiz = iz;
                    iia = ia;
                    /* Wrap around azimuth angle */
                    if (iiz < 0) {
                        iiz = nz + iiz;
                        if (dip)
                            iia = na - iia - 1;
                    } else if (iiz >= nz) {
                        iiz = iiz - nz;
                        if (dip)
                            iia = na - iia - 1;
                    }
                    gatherout[j][ia] += norm ? gatherin[iiz][iia]/(float)(2*wd + 1)
                                             : gatherin[iiz][iia];
                } /* Angles in input */
            } /* Summation base */
        }
        sf_floatwrite (gatherout[0], (nz/nth)*na, rgath);
    } /* Depths and gathers */

    free (gatherin[0]);
    free (gatherin);
    free (gatherout[0]);
    free (gatherout);

    return 0;
}

