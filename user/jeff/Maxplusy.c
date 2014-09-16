/* Computes a*x + y, where x and y are datasets, and a is scalar
x and y are floats or sf_complex, single precision
x is the stdin
a is double precision
Computations are done in double precision. */
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
    sf_file in=NULL, yfile=NULL; /* Inputs */
    sf_file out=NULL; /* Output */
    sf_file fa=NULL; /*optional input file */
    char buf_x[BUFSIZ], buf_y[BUFSIZ], buf_out[BUFSIZ]; /* I/O memory buffers */
    off_t n[SF_MAX_DIM], nsiz; /* Dims of cube, total nr elems */
    size_t i, nbuf, nleft, dim;
    size_t bufsiz;
    double a;
    float af;
    double x, y, sum; /* For the float case */
    double sum_r, sum_i, x_r, x_i, y_r, y_i; /* For the sf_complex case */
    sf_complex x_c_sp, y_c_sp; /* Single precision */
    sf_datatype type;
    bool verb;

    sf_init (argc,argv);

    in    = sf_input("in");
    yfile = sf_input("y");
    out   = sf_output("out");

    if (!sf_getdouble("a",&a)) a=1; /* Scaling factor */
    if(a==1) {
      fa = sf_input("afile");
      sf_floatread((float*)&af,1,fa   );
      a = (double)af; 
    }  






    if (!sf_getbool("verb",&verb)) verb=false; /* Verbosity flag */

    dim = (size_t) sf_largefiledims (in,n); /* Vector with cube dimensions */

    /* Total number of elements in cube as product of dimensions */
    for (nsiz=1, i=0; i < dim; i++) {
	nsiz *= n[i];
    }

    bufsiz = BUFSIZ / sf_esize(in); /* Nr of elements in buffer */
    type = sf_gettype (in);

    for (nleft=nsiz; nleft > 0; nleft -= nbuf) {
	nbuf = (bufsiz < nleft)? bufsiz: nleft;

        if (type == SF_FLOAT) {
	    sf_floatread((float*) buf_x,nbuf,in   );
            sf_floatread((float*) buf_y,nbuf,yfile);
            for (i=0; i < nbuf; i++) {
                x = ((float*)buf_x)[i];
                y = ((float*)buf_y)[i];
                sum = x * a + y;
                if(verb) sf_warning("Double precision: %1.15E", sum);
                ((float*)buf_out)[i] = (float)sum;
            }
            sf_floatwrite((float*) buf_out, nbuf, out);
        }

        else if (type == SF_COMPLEX) {
	    sf_complexread((sf_complex*) buf_x,nbuf,in);
            sf_complexread((sf_complex*) buf_y,nbuf,yfile);
            for (i=0; i < nbuf; i++) {
                x_c_sp = ((sf_complex*)buf_x)[i];
                y_c_sp = ((sf_complex*)buf_y)[i];
                x_r = (double)crealf(x_c_sp);
                y_r = (double)crealf(y_c_sp);
                x_i = (double)cimagf(x_c_sp);
                y_i = (double)cimagf(y_c_sp);
                sum_r = a * x_r + y_r;
                sum_i = a * x_i + y_i;
                if(verb) sf_warning("Double precision: %1.15E + %1.15E i", sum_r, sum_i);
                ((float*)buf_out)[2*i  ] = (float)sum_r;
                ((float*)buf_out)[2*i+1] = (float)sum_i;
            }
            sf_floatwrite((float*) buf_out, 2*nbuf, out);
	}
    }

    exit (0);
}
