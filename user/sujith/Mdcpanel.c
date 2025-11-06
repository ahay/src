/* Phase velocity vs frequency panels for surface waves */
/*
  Copyright (C) 2025 University of Texas at Austin

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
#include <complex.h>

int main(int argc, char* argv[])
{
    int nf, nx, n3, if0, ix, nc, i3, ic;
    float of, df, ox, dx, oc, dc;
    sf_file in, out;

    float *buf = NULL;     
    float *panel = NULL;  
    float *x = NULL;
    float *c = NULL;
    float *f = NULL;

    sf_init (argc, argv);

    in  = sf_input ("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input (SF_COMPLEX)");

    if (!sf_histint(in,"n1",&nf)) sf_error("No n1=");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2=");
    if (!sf_histfloat(in,"o1",&of)) sf_error("No o1=");
    if (!sf_histfloat(in,"d1",&df)) sf_error("No d1=");
    if (!sf_histfloat(in,"o2",&ox)) sf_error("No o2=");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2=");

    n3 = sf_leftsize(in,2);

    
    if (!sf_getfloat("c0",&oc)) oc = 500.0f;    /* starting phase velocity */
    if (!sf_getfloat("dc",&dc)) dc = 10.0f;     /* phase velocity increment */
    if (!sf_getint("nc",&nc)) nc = 100;         /* number of phase-velocity samples */

    /* allocate buffers */
    buf   = sf_floatalloc(2 * nf * nx);   
    panel = sf_floatalloc(nf * nc);       
    x     = sf_floatalloc(nx);
    c     = sf_floatalloc(nc);
    f     = sf_floatalloc(nf);

    /* write out header for a real output */
    sf_putfloat(out,"o1",of);
    sf_putfloat(out,"d1",df);
    sf_putint(out,"n1",nf);
    sf_putfloat(out,"o2", oc);
    sf_putfloat(out,"d2", dc);
    sf_putint(out,"n2", nc);
    sf_settype(out, SF_FLOAT);

    /* coords */
    for (ix = 0; ix < nx; ix++) x[ix] = ox + ix * dx;
    for (ic = 0; ic < nc; ic++) c[ic] = oc + ic * dc;
    for (if0 = 0; if0 < nf; if0++) f[if0] = of + if0 * df;

    sf_warning("Min freq: %.2f, Max freq: %.2f", f[0], f[nf-1]);
    sf_warning("Min offset: %.2f, Max offset: %.2f", x[0], x[nx-1]);
    sf_warning("Min phase vel: %.2f, Max phase vel: %.2f", c[0], c[nc-1]);

    /* main loop over 3rd dimension */
    for (i3 = 0; i3 < n3; i3++) {

        sf_floatread(buf, 2 * nf * nx, in);

        /* compute panel */
        for (if0 = 0; if0 < nf; if0++) {
            for (ic = 0; ic < nc; ic++) {
                float k = 2.0f * M_PI * f[if0] / c[ic]; /* w / c */

                float complex sum = 0.0f + 0.0f * I;

                /* sum over space */
                for (ix = 0; ix < nx; ix++) {

                    int base = (if0 + nf * ix) * 2;
                    float real_part = buf[base];
                    float imag_part = buf[base + 1];
                    float magnitude = sqrtf(real_part * real_part + imag_part * imag_part) + FLT_EPSILON;
                    real_part /= magnitude;
                    imag_part /= magnitude;

                    float complex samp = real_part + imag_part * I;

                    float complex phase = cexpf(I * k * x[ix]);
                    sum += samp * phase;
                }

                panel[if0 + nf * ic] = cabsf(sum);
            }
        }

        /* write panel to output */
        sf_floatwrite(panel, nf * nc, out);
    }

    free(buf);
    free(panel);
    free(x);
    free(c);
    free(f);

    return 0;
}
