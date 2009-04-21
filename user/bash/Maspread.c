/* Spread surface(s) with analytic expressions in a data volume

Example:
sfaspread n1=501 n2=101 n3=101 o2=-0.5 o3=-0.5 d2=0.01 d3=0.01
          nsp=2 pos1="sqrt (0.25 + x2*x2 + x3*x3)"
                pos2="1.5 + 0.5*x2 + 0.5*x3" > out.rsf

nsp defines number of surfaces, pos# defines an expression
for the #th surface, amp# defines an expression for the
amplitude of the #th surface (it is optional and assumed
to be "1" by default).

See sfmath for the details about the syntax of expressions. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

static float sf_ricker_impulse (float f0, float tau, float t) {
    return (1.0 - 2.0 * SF_PI * SF_PI * f0 * f0 * (tau - t) * (tau - t)) *
           exp(-SF_PI * SF_PI * f0 * f0 * (tau - t) * (tau - t));
}

int main (int argc, char* argv[]) {
    const unsigned int max_buf = 1 << 24;
    char* unit_amplitude = "1";
    int i, k, l, dim, n[SF_MAX_DIM], ii[SF_MAX_DIM], nbuf;
    off_t j, nsiz;
    size_t len;
    int nsp, n1, t0, t1;
    float **fbuf, ***p, ***a, d[SF_MAX_DIM], o[SF_MAX_DIM], *trace;
    float f, fpeak, d1, o1;
    char key[7], *label, *unit, **pos, **amp;
    bool conv, grad;
    sf_file out;

    sf_init (argc, argv);
    out = sf_output ("out");
    sf_setformat (out, "native_float");

    if (!sf_getbool ("conv", &conv)) conv = true;
    /* y - convolve every surface with Ricker wavelet */

    if (!sf_getbool ("grad", &grad)) grad = false;
    /* y - use gradient to fill the space below every
           surface with amplitude value (if conv=y, it is ignored) */

    /* dimensions */
    for (i = 0; i < SF_MAX_DIM; i++) {
        snprintf (key, 3, "n%d", i + 1);
        if (!sf_getint (key, n + i)) break;
        /*( n# size of #-th axis )*/
        sf_putint (out, key, n[i]);
    }

    if (0 == i) sf_error ("Need n1=");
    dim = i;
    n1 = n[0];

    /* basic parameters */
    for (i = 0; i < dim; i++) {
        snprintf (key, 3, "o%d", i + 1);
        if (!sf_getfloat (key, &f)) f = 0.;
        /*( o#=(0,0,...) origin on #-th axis )*/
        sf_putfloat (out, key, f);
        o[i] = f;

        snprintf (key, 3, "d%d", i + 1);
        if (!sf_getfloat (key, &f)) f = (i == 0) ? 0.004 : 0.1;
        /*( d#=(0.004,0.1,0.1,...) sampling on #-th axis )*/
        sf_putfloat (out, key, f);
        d[i] = f;

        snprintf (key, 7, "label%d", i + 1);
        if (NULL == (label = sf_getstring (key)))
            label = (i == 0) ? "Time" : "Distance";
        /*( label#=(Time,Distance,Distance,...) label on #-th axis )*/
        if (*label != '\0' && (*label != ' ' || *(label + 1) != '\0'))
            sf_putstring (out, key, label);

        snprintf (key, 6, "unit%d", i + 1);
        if (NULL == (unit = sf_getstring (key)))
            unit = (i == 0) ? "s" : "km";
        /*( unit#=[s,km,km,...] unit on #-th axis )*/
        if (*unit != '\0' && (*unit != ' ' || *(unit + 1) != '\0'))
            sf_putstring (out, key, unit);
    }
    o1 = o[0];
    d1 = d[0];

    if (NULL != (label = sf_getstring ("title")))
        sf_putstring (out, "title", label);
    /* title for plots */

    if (!sf_getfloat ("fpeak", &fpeak)) fpeak = 20;
    /* peak frequency for Ricker wavelet (in Hz) */

    if (!sf_getint ("nsp", &nsp) || nsp < 1) nsp = 1;
    /* Number of surfaces */

    pos = (char**)malloc (nsp * sizeof (char*));
    amp = (char**)malloc ((nsp + 2) * sizeof (char*));

    /* expressions for positions and amplitudes */
    amp[0] = NULL;
    for (i = 0; i < nsp; i++) {
        pos[i] = NULL;
        amp[i + 1] = NULL;
        snprintf (key, 5, "pos%d", i + 1);
        if (NULL == (pos[i] = sf_getstring (key))) sf_error ("Need %s=", key);
        /* Describes the vertical position in a mathematical notation. */
        snprintf (key, 5, "amp%d", i + 1);
        if (NULL == (amp[i + 1] = sf_getstring (key))) amp[i + 1] = unit_amplitude;
        /* Describes the amplitude in a mathematical notation. */
    }
    if (!conv) {
       snprintf (key, 5, "amp%d", 0);
       if (NULL == (amp[0] = sf_getstring (key))) amp[0] = unit_amplitude;
       /* If conv=false, then this describes values at the top. */
       if (grad) {
          snprintf (key, 5, "amp%d", nsp + 1);
          if (NULL == (amp[nsp + 1] = sf_getstring (key))) amp[nsp + 1] = amp[nsp];
          /* If grad=true, then this describes values at the bottom. */
       }
    }

    for (i = 0; i < dim; i++) {
	sprintf (key, "x%d", i + 1);

	if (NULL != sf_histstring (out, key))
	    sf_error ("illegal use of %s parameter", key);
	sf_putint (out, key, i);

        ii[i] = 0;
    }

    trace = sf_floatalloc (n1);
    /* Fale write to flush the header */
    sf_floatwrite (trace, 0, out);
    n[0] = 1;
    o[0] = 1;
    d[0] = 1;
    /* Put fake values for math parser */
    sf_putint (out, "n1", 1);
    sf_putfloat (out, "o1", 1);
    sf_putfloat (out, "d1", 1);

    p = (float***)malloc (nsp * sizeof (float**));
    a = (float***)malloc ((nsp + 2) * sizeof (float**));

    nbuf = max_buf / (sizeof (float) * 2 * nsp * dim);
    fbuf = sf_floatalloc2 (nbuf, dim);

    for (nsiz = 1, i = 0; i < dim; i++) {
	nsiz *= n[i];
    }

    /* Loop over buffer chunks */
    for (j = 0; nsiz > 0; nsiz -= nbuf) {
	if (nbuf > nsiz) nbuf = nsiz;

        /* Fill out coordinates */
	for (k = 0; k < nbuf; k++, j++) {
	    sf_line2cart (dim, n, j, ii);
	    for (i = 0; i < dim; i++) {
		fbuf[i][k] = o[i] + ii[i] * d[i];
	    }
	}

        /* Loop over surfaces */
        for (k = 0; k < nsp; k++) {
            len = sf_math_parse (pos[k], out, SF_FLOAT);
            p[k]  = sf_floatalloc2 (nbuf, len + 2);
	    sf_math_evaluate (len, nbuf, fbuf, p[k]);

            len = sf_math_parse (amp[k + 1], out, SF_FLOAT);
            a[k + 1]  = sf_floatalloc2 (nbuf, len + 2);
	    sf_math_evaluate (len, nbuf, fbuf, a[k + 1]);
        }
        if (!conv) {
            len = sf_math_parse (amp[0], out, SF_FLOAT);
            a[0]  = sf_floatalloc2 (nbuf, len + 2);
	    sf_math_evaluate (len, nbuf, fbuf, a[0]);
            if (grad) {
                len = sf_math_parse (amp[nsp + 1], out, SF_FLOAT);
                a[nsp + 1]  = sf_floatalloc2 (nbuf, len + 2);
	        sf_math_evaluate (len, nbuf, fbuf, a[nsp + 1]);
            }
        }

        if (conv) { /* Convolve with Ricker */
            /* Loop over traces */
	    for (k = 0; k < nbuf; k++) {
                /* Loop over time steps */
                for (i = 0; i < n1; i++) {
                     trace[i] = 0.0;
                     /* Loop over surfaces */
                     for (l = 0; l < nsp; l++)
                         trace[i] += a[l + 1][1][k] *
                                     sf_ricker_impulse (fpeak, o1 + i * d1,
                                                               p[l][1][k]);
                }
	        sf_floatwrite (trace, n1, out);
            }
        } else { /* Fill values below each surface */
            /* Loop over traces */
	    for (k = 0; k < nbuf; k++) {
                /* Loop over surfaces */
                t1 = 0;
                for (l = 0; l <= nsp; l++) {
                    f = 0;
                    t0 = t1;
                    t1 = l < nsp ? (int)((p[l][1][k] - o1) / d1 + 0.5) : n1;
                    if (t1 < 0)
                        t1 = 0;
                    if (t1 > n1)
                        t1 = n1;
                    if (grad && t0 != t1)
                        f = (a[l + 1][1][k] - a[l][1][k]) / (t1 - t0);
                    /* Loop over time steps */
                    for (i = t0; i != t1; t1 > t0 ? i++ : i--) {
                        trace[i] = a[l][1][k] + f * (i - t0);
                    }
                }
	        sf_floatwrite (trace, n1, out);
            }
        }
        for (k = 0; k < nsp; k++) {
            free (p[k]);
            free (a[k + 1]);
        }
        if (!conv) {
            free (a[0]);
            if (grad)
                free (a[nsp + 1]);
        }
    }

    exit (0);
}
