/* Compute RBF across fault using fault attribute computed by Sobel filter. */
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
#include <math.h>
#include <stdlib.h>

void normalize(int n1, int n2, float **scan, float scaleFac);

int main(int argc, char *argv[])
{
    int n1,n2,n3,n12; // n1 is trace time length, n2 is trace number
    int i1,i2,i3,j2; // index on each dimension
    int iRef; // [input] reference trace index
    float **sobel; // [input] n1*n2 Sobel filter fault attribute
    float **rbf; // [output] n1*n2 computed RBF across fault
    float scaleFac; // [input] scale factor
    float rRef; // [input] reference radial in RBF
    float radial; // radial varialbe used in RBF
    bool useinput; // [input] bool flag: whether use input fault attribute

    sf_file in,out;
    sf_init(argc,argv);
    in = sf_input("in"); // input Sobel fault attribute
    /* DEFAULT: input */
    out = sf_output("out"); // output RBF across fault
    /* DEFAULT: output */
    if (!sf_getint("i0",&iRef)) sf_error("Need reference trace i0.");
    /* Reference trace position. */
    if (!sf_getfloat("scale",&scaleFac)) scaleFac=1.0;
    /* Fault attribute scaling factor (0.0 ~ factor). */
    if (!sf_getfloat("r0",&rRef)) rRef=1.0;
    /* Reference radial in RBF. */
    if (!sf_getbool("useinput",&useinput)) useinput=true;
    /* Flag: whether use the input fault attribute. */

    // Check input data dimensions
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input.");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input.");
    n12 = n1 * n2;
    n3 = sf_leftsize(in,2);
    if (iRef<0 || iRef >= n2) sf_error("Reference i0 out of boundary.");

    // Memory allocation for all arrays
    sobel = sf_floatalloc2(n1,n2);
    rbf = sf_floatalloc2(n1,n2);

    // Most outside loop on n3
    for (i3=0; i3<n3; i3++) {
        if (useinput) {
            // Read in input file
            sf_floatread(sobel[0],n12,in);
            // Apply fault attribute normalization
            normalize(n1,n2,sobel,scaleFac);
        } else {
            for (i2=0; i2<n2; i2++)
                for (i1=0; i1<n1; i1++)
                    sobel[i2][i1] = 0.0;
        }

        // Initialize other arrays
        for (i2=0; i2<n2; i2++)
            for (i1=0; i1<n1; i1++)
                rbf[i2][i1] = 0.0;

        // Essential computations
        for (i2=0; i2<n2; i2++) {
            for (i1=0; i1<n1; i1++) {
                radial=0.0;
                if (i2 < iRef) {
                    for (j2=i2+1; j2<=iRef; j2++) radial += exp(sobel[j2][i1]);
                    // for (j2=i2+1; j2<=iRef; j2++) radial += 1.0;
                } else if (i2 == iRef) {
                    radial = 0.0;
                } else { // if (i2 > iRef)
                    for (j2=iRef; j2<i2; j2++) radial += exp(sobel[j2][i1]);
                    // for (j2=iRef; j2<i2; j2++) radial += 1.0;
                }

                radial /= rRef;
                rbf[i2][i1] = 1.0 / (1.0 + radial*radial);
            }
        }

        // Write into output files
        sf_floatwrite(rbf[0],n12,out);
    }

    exit(0);
}

void normalize(int n1, int n2, float **scan, float scaleFac) {
    int i1,i2;
    float min=FLT_MAX,max=0.0;

    for (i1=0; i1 < n1; i1++) {
        for (i2=0; i2 < n2; i2++) {
            if (scan[i2][i1] < min) min = scan[i2][i1];
            if (scan[i2][i1] > max) max = scan[i2][i1];
        }
    }
    if ((max-min) < 1e-6) sf_warning("WARNING: Input fault range < 1e-6.");

    for (i1=0; i1 < n1; i1++)
        for (i2=0; i2 < n2; i2++)
            scan[i2][i1] = (scan[i2][i1]-min) / (max-min) * scaleFac;
}
