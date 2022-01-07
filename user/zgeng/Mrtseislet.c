/* Seislet transform using relative time */
/*
  Copyright (C) 2018 University of Texas at Austin
   
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

#include "rtseislet.h"

int main(int argc, char *argv[])
{
    int n1, n2, n3, n12, nt, nw; // n1: time length, n2: trace number
    float o1, d1;
    int i1, i2, i3; // auxiliary index
    char *type;     // wavelet type
    int order;
    bool inv, adj, unit, verb;
    float *pp, *qq, **str, eps, **sinv, **time;
    sf_map4 mo;
    sf_file in, out, rt;

    sf_init(argc, argv);

    in  = sf_input("in");   // input traces file
    out = sf_output("out"); // output seislet file
    rt  = sf_input("rt");   // input rt file

    if (!sf_histint(in, "n1", &n1))
        sf_error("No n1= in input");
    if (!sf_histint(in, "n2", &n2))
        sf_error("No n2= in input");
    n3 = sf_leftsize(in, 2);
    n12 = n1 * n2;

    if (!sf_histfloat(in, "d1", &d1))
        d1 = 1.;
    if (!sf_histfloat(in, "o1", &o1))
        o1 = 0.;

    if (!sf_histint(rt, "n1", &nt))
        sf_error("No n1= in rt");
    if (nt != n1)
        sf_error("n1=%d in rt must equal n1=%d in input", nt, n1);
    if (!sf_histint(rt, "n2", &nw))
        sf_error("No n2= in rt");
    if (nw != n2)
        sf_error("Need %d traces in rt, got %d", n2, nw);

    pp = sf_floatalloc(n12);
    qq = sf_floatalloc(n12);
    str = sf_floatalloc2(n1, n2);
    sinv = sf_floatalloc2(n1, n2);
    time = sf_floatalloc2(n1, n2);

    if (!sf_getbool("inv", &inv))
        inv = false;
    /* if y, do inverse transform */

    if (!sf_getbool("adj", &adj))
        adj = false;
    /* if y, do adjoint transform */

    if (!sf_getfloat("eps", &eps))
        eps = 0.01;
    /* regularization */

    if (!sf_getbool("unit", &unit))
        unit = false;
    /* if y, use unitary scaling */

    if (!sf_getint("order", &order))
        order = 1;
    /* accuracy order */

    if (NULL == (type = sf_getstring("type")))
        type = "linear";
    /* [haar,linear,biorthogonal] wavelet type, the default is linear  */

    if (!sf_getbool("verb", &verb))
        verb = false;
    /* verbosity flag */

    /* Initialize for seislet transform */
    seislet_init(n1, n2, o1, d1, inv, unit, eps, order, type[0]);
    seislet_set(str, sinv);

    mo = sf_stretch4_init(n1, o1, d1, n1, 1);

    for (i3 = 0; i3 < n3; i3++)
    {
        if (verb)
            sf_warning("slice %d of %d;", i3 + 1, n3);

        sf_floatread(pp, n12, in);
        sf_floatread(str[0], n1 * n2, rt);

        /* Inverse interpolate T0(x,tau)*/
        for (i2 = 0; i2 < n2; i2++)
        {
            for (i1 = 0; i1 < n1; i1++)
                time[i2][i1] = o1 + d1 * i1;
            sf_stretch4_define(mo, str[i2],false);
            sf_stretch4_apply(false, mo, time[i2], sinv[i2]);
        }

        /* Clean shift */
        for (i2 = 0; i2 < n2; i2++)
            for (i1 = 0; i1 < n1 - 1; i1++)
            {
                if (str[i2][i1 + 1] - str[i2][i1] < 0.001)
                {
                    str[i2][i1 + 1] = str[i2][i1] + 0.001;
                }
            }

        /* Seislet transform */
        if (adj)
        {
            seislet_lop(adj, false, n12, n12, qq, pp);
        }
        else
        {
            seislet_lop(adj, false, n12, n12, pp, qq);
        }

        /* Write into output file */
        sf_floatwrite(qq, n12, out);
    }
    exit(0);
}
