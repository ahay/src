/* Transforms regular 2-D array to sparse array
Input is int or float array
Output is a float nonzero-by-3 array, where 
nonzero=`<input.rsf sfattr want=nonzero | awk -F '= ' '{ print $2 }';`
column 0 holds the data values (converted from int to float, if necessary),
column 1 holds coordinate values (i.e. o+i*d, not indices i) for dimension
1 of input, and column 2 the coordinate values for dimension 2 of input */
/*
  Copyright (C) 2011 Ioan Vlad

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

int main(int argc, char* argv[]) 
{
    int i1, i2;  /* Counters over input */
    int n1, n2, nonzero;  /* Input dims, nr of nonzero elements */
    int i; /* Counter over output */
    int n; /* Output size */
    int ndims_in = 2; /* Input dimensionality */
    int n2_out; /* Output n2 */
    int   *trci=NULL; /* Input trace (int data) */
    float *trcf=NULL; /* Input trace (float data) */
    float *oa=NULL; /* Output array */
    float o1, d1, o2, d2; /* Input axes origin and sampling */
    sf_file in=NULL, out=NULL; /* Input, output files */
    sf_datatype inp_type; /* Input data type */

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in, "n1", &n1)) sf_error("Need n1=");
    /* number of histogram samples */
    if (!sf_histfloat(in, "o1", &o1)) sf_error("Need o1=");
    /* histogram origin */
    if (!sf_histfloat(in, "d1", &d1)) sf_error("Need d1=");
    /* histogram sampling */
    if (!sf_histint(in, "n2", &n2)) sf_error("Need n2=");
    /* number of histogram samples in dimension 2 */
    if (!sf_histfloat(in, "o2", &o2)) sf_error("Need o2=");
    /* histogram origin for dimension 2 */
    if (!sf_histfloat(in, "d2", &d2)) sf_error("Need d2=");
    /* histogram sampling for dimension 2 */

    if (!sf_getint("nonzero",&nonzero)) nonzero = n1*n2;
    /* Number of nonzero elements in input */

    inp_type = sf_gettype(in);

    if (inp_type == SF_INT) {
        trci = sf_intalloc(n1);
        sf_settype(out,SF_FLOAT);
    }
    else if (inp_type != SF_FLOAT) {
        sf_error("Need float or int input");
    }

    n2_out = ndims_in + 1;
    n = n2_out * nonzero;
    trcf = sf_floatalloc(n1);
    oa = sf_floatalloc(n);

    sf_putint(out, "n1", nonzero);
    sf_putint(out, "n2", n2_out);

    i = 0;

    /* Duplicating boilerplate code to avoid conditionals inside loops */

    for (i2=0; i2 < n2; i2++) {

        if (inp_type == SF_INT) {
            sf_intread(trci, n1, in);
            for (i1=0; i1 < n1; i1++) {
                trcf[i1] = trci[i1];
            }
        } else {
            sf_floatread(trcf,n1,in);
        }

        for (i1=0; i1 < n1; i1++) {
            if (trcf[i1] != 0) {
                if (i==nonzero) sf_error("nonzero was too small!");
                oa[i            ] = trcf[i1];
                oa[i + nonzero  ] = o1 + i1*d1;
                oa[i + nonzero*2] = o2 + i2*d2;
                i++;
            }
        }
    }
    if (i<nonzero-1) sf_error("nonzero was too big!");
    sf_floatwrite(oa, n, out);

    exit(0);
}
