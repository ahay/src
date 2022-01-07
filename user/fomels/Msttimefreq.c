/* Time-frequency analysis using streaming attributes. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include "stdiv.h"

int main(int argc, char *argv[])
{
    bool phase;
    char *label;
    int i1, n1, iw, nt, nw, i2, n2, window_size;
    float t, d1, w, w0, dw, *trace, *bs, *bc, *ss, *cc, *mm, lam;
    sf_file time, timefreq, mask;

    sf_init(argc, argv);
    time = sf_input("in");
    timefreq = sf_output("out");

    if (!sf_histint(time, "n1", &n1))
        sf_error("No n1= in input");
    if (!sf_histfloat(time, "d1", &d1))
        d1 = 1.;
    n2 = sf_leftsize(time, 1);

    if (!sf_getint("nw", &nw))
    { /* number of frequencies */
        nt = 2 * kiss_fft_next_fast_size((n1 + 1) / 2);
        nw = nt / 2 + 1;
        dw = 1. / (nt * d1);
        w0 = 0.;
    }
    else
    {
        if (!sf_getfloat("dw", &dw))
        {
            /* f	requency step */
            nt = 2 * kiss_fft_next_fast_size((n1 + 1) / 2);
            dw = 1. / (nt * d1);
        }
        if (!sf_getfloat("w0", &w0))
            w0 = 0.;
        /* first frequency */
    }

    sf_shiftdim(time, timefreq, 2);

    sf_putint(timefreq, "n2", nw);
    sf_putfloat(timefreq, "d2", dw);
    sf_putfloat(timefreq, "o2", w0);

    if (NULL != (label = sf_histstring(time, "label1")) && !sf_fft_label(2, label, timefreq))
        sf_putstring(timefreq, "label2", "Wavenumber");
    sf_fft_unit(2, sf_histstring(time, "unit1"), timefreq);

    dw *= 2. * SF_PI;
    w0 *= 2. * SF_PI;

    trace = sf_floatalloc(n1);
    bs = sf_floatalloc(n1);
    bc = sf_floatalloc(n1);
    ss = sf_floatalloc(n1);
    cc = sf_floatalloc(n1);

    if (!sf_getfloat("lambda", &lam))
        lam = 1.0f;
    /* smoothing parameter */

    if (!sf_getint("window_size", &window_size))
        window_size = 50;
    /* window_size */

    stdiv_init(n1, window_size, lam);

    if (!sf_getbool("phase", &phase))
        phase = false;
    /* output phase instead of amplitude */

    if (NULL != sf_getstring("mask"))
    {
        mask = sf_input("mask");
        mm = sf_floatalloc(n1);
    }
    else
    {
        mask = NULL;
        mm = NULL;
    }

    for (i2 = 0; i2 < n2; i2++)
    {
        sf_floatread(trace, n1, time);
        if (NULL != mm)
        {
            sf_floatread(mm, n1, mask);
            for (i1 = 0; i1 < n1; i1++)
            {
                trace[i1] *= mm[i1];
            }
        }

        for (iw = 0; iw < nw; iw++)
        {
            w = w0 + iw * dw;

            if (0. == w)
            { /* zero frequency */
                for (i1 = 0; i1 < n1; i1++)
                {
                    ss[i1] = 0.;
                    bc[i1] = 0.5;
                    if (NULL != mm)
                        bc[i1] *= mm[i1];
                }

                stdiv(trace, bc, cc);
            }
            else
            {
                for (i1 = 0; i1 < n1; i1++)
                {
                    t = i1 * d1;
                    bs[i1] = sinf(w * t);
                    bc[i1] = cosf(w * t);
                    if (NULL != mm)
                    {
                        bs[i1] *= mm[i1];
                        bc[i1] *= mm[i1];
                    }
                }

                stdiv(trace, bs, ss);
                stdiv(trace, bc, cc);
            }

            for (i1 = 0; i1 < n1; i1++)
            {
                ss[i1] = phase ? atan2f(ss[i1], cc[i1]) : hypotf(ss[i1], cc[i1]);
            }

            sf_floatwrite(ss, n1, timefreq);
        }
    }

    exit(0);
}
