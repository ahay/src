/* Mask values inside or outside of a range.

sfmaskval < in.rsf > out.rsf upper=ul lower=ll upperval=uv lowerval=lv

If upper > lower, then values larger than ul will be changed to uv and
values belowe ll will be changed to lv.
If upper < lower, then values inside [ul;ll] will be changed to lv.

*/

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
/* Based upon sfclip2 */


#include <rsf.h>
#include <float.h>

int main (int argc, char* argv[]) {
    int i, n, nbuf;
    float upper, lower, upperval, lowerval, *trace=NULL;
    sf_file in, out;

    /* Initialize RSF */
    sf_init (argc, argv);
    /* standard input */
    in = sf_input ("in");
    /* standard output */
    out = sf_output ("out");

    /* check that the input is float */
    if (SF_FLOAT != sf_gettype (in)) sf_error ("Need float input");

    n = sf_filesize (in);

    if (!sf_getfloat ("upper", &upper)) upper = +FLT_MAX;
    /* upper range limit */
    if (!sf_getfloat ("lower", &lower)) lower = -FLT_MAX;
    /* lower range limit */
    if (!sf_getfloat ("upperval", &upperval)) upperval = +FLT_MAX;
    /* upper range value */
    if (!sf_getfloat ("lowerval", &lowerval)) lowerval = -FLT_MAX;
    /* lower range value */

    nbuf = BUFSIZ/sizeof(float);
    trace = sf_floatalloc (nbuf);

    /* loop over traces */
    for (n = sf_filesize (in); n > 0; n -= nbuf) {
        if (nbuf > n) nbuf=n;

        sf_floatread (trace, nbuf, in);

        for (i = 0; i < nbuf; i++) {
            if (upper > lower) {
                if (trace[i] > upper)
                    trace[i] = upperval;
                if (trace[i] < lower)
                    trace[i] = lowerval;
            } else {
                if (trace[i] > upper &&
                    trace[i] < lower)
                    trace[i] = lowerval;
            }
        }

        sf_floatwrite (trace, nbuf, out);
    }

    return 0;
}

