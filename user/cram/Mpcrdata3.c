/* Prepare data for 3-D angle-domain migration. */
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
    int i, n, nt;
    float *tracein, *traceout, dt;
    sf_file data, out;

    bool filter, kmah, diff, erefl, verb;

    sf_init (argc, argv);

    data = sf_input ("in");
    /* Common-shot 2-D data */
    out = sf_output ("out");
    /* Traces for migration */

    if (!sf_getbool ("verb", &verb)) verb = false;
    /* verbosity flag */

    if (!sf_getbool ("erefl", &erefl)) erefl = false;
    /* y - assume data modeled with exploding reflector */
    if (!sf_getbool ("filter", &filter)) filter = true;
    /* y - antialiasing filter for data */
    if (!sf_getbool ("KMAH", &kmah)) kmah = true;
    /* y - account for phase shifts due to KMAH index */
    if (!sf_getbool ("diff", &diff)) diff = true;
    /* y - apply differentiation */

    /* Data dimensions */
    if (!sf_histint (data, "n1", &nt)) sf_error ("No n1= in data");
    if (!sf_histfloat (data, "d1", &dt)) sf_error ("No d1= in data");
    n = sf_leftsize (data, 1);

    tracein = sf_floatalloc (nt);
    traceout = sf_floatalloc (nt*2);

    if (kmah) { /* Make room for phase-shifted traces */
        sf_shiftdim (data, out, 2);
        sf_putint (out, "n2", 2);
        sf_putfloat (out, "o2", 0.0);
        sf_putfloat (out, "d2", 1.0);
        sf_putstring (out, "label1", "");
        sf_putstring (out, "unit1", "");
    }
    sf_putstring (out, "ExplRefl", erefl ? "y" : "n");
    sf_putstring (out, "filter", filter ? "y" : "n");
    sf_putstring (out, "KMAH", kmah ? "y" : "n");

    for (i = 0; i < n; i++) { /* Loop over traces */
        if (verb && !(i % 10000LU))
            sf_warning ("Processing trace %lu of %lu (%g %%);",
                        i + 1LU, n, 100.0*(float)i/(float)n);
        sf_floatread (diff ? tracein : traceout, nt, data);
        /* Differentiate */
        if (diff)
            sf_cram_trace_deriv (tracein, traceout, nt, dt);
        /* pi/2 phase shift */
        if (kmah)
            sf_cram_trace_hilbert (nt, traceout, &traceout[nt]);
        /* Causal and anti-causal integration for anti-aliasing filter */
        if (filter) {
            sf_cram_trace_cint (traceout, nt);
            sf_cram_trace_acint (traceout, nt);
            if (kmah) {
                sf_cram_trace_cint (&traceout[nt], nt);
                sf_cram_trace_acint (&traceout[nt], nt);
            }
        }
        sf_floatwrite (traceout, nt, out);
        if (kmah)
            sf_floatwrite (&traceout[nt], nt, out);
    }
    if (verb)
        sf_warning (".");

    free (tracein);
    free (traceout);

    return 0;
}

