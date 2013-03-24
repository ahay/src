/* Computes square of L-2 norm in double precision.
*/
/*
  Copyright (C) 2004 University of Texas at Austin
  Copyright (C) 2009 Colorado School of Mines
  Copyright (C) 2009 Statoil ASA

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

#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <rsf.h>

int main(int argc, char* argv[])
{
    sf_file in=NULL;
    char buf[BUFSIZ];
    off_t n[SF_MAX_DIM], nsiz;
    size_t i, nbuf, nleft, dim;
    size_t bufsiz;
    float f;
    double f_double, c_real, c_imag;
    double fsqr;
    sf_complex c;
    sf_datatype type;

    sf_init (argc,argv);
    in = sf_input("in");

    dim = (size_t) sf_largefiledims (in,n); /* Vector with cube dimensions */

    /* Total number of elements in cube as product of dimensions */
    for (nsiz=1, i=0; i < dim; i++) {
	nsiz *= n[i];
    }

    bufsiz = BUFSIZ / sf_esize(in); /* Nr of elements in buffer */
    type = sf_gettype (in);
    fsqr = 0; /* Summation result */

    /* A little bit of code duplication in order to avoid testing the type
    for every single value in the dataset. The clean way to do this would be
    to have a function that takes another function as an argument */

    switch (type) {
        case SF_FLOAT:
            for (nleft=nsiz; nleft > 0; nleft -= nbuf) {
                nbuf = (bufsiz < nleft)? bufsiz: nleft;
                /* Read nbuf elements from "in" into buf. */
                sf_floatread((float*) buf,   nbuf, in);
                for (i=0; i < nbuf; i++) {
                    f        = ((float*)buf)[i];
                    f_double = (double) f;
                    fsqr += f_double*f_double;
                }
            }
        break;
        case SF_COMPLEX:
            for (nleft=nsiz; nleft > 0; nleft -= nbuf) {
                nbuf = (bufsiz < nleft)? bufsiz: nleft;
                /* Read nbuf elements from "in" into buf. */
                sf_complexread((sf_complex*) buf,  nbuf,in);
                for (i=0; i < nbuf; i++) {
                    c = ((sf_complex*)buf)[i];
                    c_real = (double) crealf(c);
                    c_imag = (double) cimagf(c);
                    fsqr += c_real*c_real + c_imag*c_imag;
                }
            }
        break;
	default:
	    sf_error("Bad type %d",type);
	    break;
    }

    printf("%1.15E\n", fsqr);

    exit (0);
}
