/* Prepare data for 2-D angle-domain migration. */
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

#include "cram_helper.h"

int main (int argc, char* argv[]) {
    int i, j, n, nn, nt;
    float *trace;
    sf_file data, out;

    bool absoff, filter, kmah, diff, verb;

    sf_init (argc, argv);

    data = sf_input ("in");
    /* Common-shot 2-D data */
    out = sf_output ("out");

    if (!sf_getbool ("verb", &verb)) verb = false;
    /* verbosity flag */

    if (!sf_getbool ("absoff", &absoff)) absoff = false;
    /* y - absolute offset (default - relative to shot axis) */
    if (!sf_getbool ("filter", &filter)) filter = true;
    /* y - antialiasing filter for data */
    if (!sf_getbool ("KMAH", &kmah)) kmah = true;
    /* y - account for phase shifts due to KMAH index */
    if (!sf_getbool ("diff", &diff)) diff = true;
    /* y - apply half-order differentiation */

    /* Data dimensions */
    if (!sf_histint (data, "n1", &nt)) sf_error ("No n1= in data");
    n = sf_leftsize (data, 1);

    /* Next suitable size for the differentiator */
    nn = 2*kiss_fft_next_fast_size ((nt + 1)/2);
    trace = sf_floatalloc (nn*2);

    if (kmah) { /* Make room for phase-shifted traces */
        sf_shiftdim (data, out, 2);
        sf_putint (out, "n2", 2);
        sf_putfloat (out, "o2", 0.0);
        sf_putfloat (out, "d2", 1.0);
        sf_putstring (out, "label1", "");
        sf_putstring (out, "unit1", "");
    }
    sf_putstring (out, "absoff", absoff ? "y" : "n");
    sf_putstring (out, "filter", filter ? "y" : "n");
    sf_putstring (out, "KMAH", kmah ? "y" : "n");

    /* Half-order differentiation object */
    if (diff)
        sf_halfint_init (true, nn, 1.-1./nt);

    for (i = 0; i < n; i++) { /* Loop over traces */
        if (verb && !(i % 10000LU))
            sf_warning ("Processing trace %lu of %lu (%g %%);",
                        i + 1LU, n, 100.0*(float)i/(float)n);
        sf_floatread (trace, nt, data);
        for (j = nt; j < nn; j++) {
            trace[j] = 0.;
        }
        /* Differentiate */
        if (diff) 
            sf_halfint (true, trace);
        /* pi/2 phase shift */
        if (kmah)
            sf_cram_trace_hilbert (nt, trace, &trace[nn]);
        /* Causal and anti-causal integration for anti-aliasing filter */
        if (filter) {
            sf_cram_trace_cint (trace, nt);
            sf_cram_trace_acint (trace, nt);
            if (kmah) {
                sf_cram_trace_cint (&trace[nn], nt);
                sf_cram_trace_acint (&trace[nn], nt);
            }
        }
        sf_floatwrite (trace, nt, out);
        if (kmah)
            sf_floatwrite (&trace[nn], nt, out);
    }
    if (verb)
        sf_warning (".");

    free (trace);

    return 0;
}

