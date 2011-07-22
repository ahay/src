/* Check whether all values in a dataset are a multiple of a factor or not */
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
#include <math.h>

int main(int argc, char* argv[]) 
{
    int i, n, nbuf;
    int glob_i = 0;   /* Global counter in file */
    int   *ibuf=NULL; /* Input reading buffer */
    float *fbuf=NULL; /* Input reading buffer */
    sf_file in=NULL;
    bool float_inp=false, int_inp=false; /* Input flags */
    sf_datatype inp_type; /* Input data type */
    int i_fac;     /* Integer division factor */
    float f_fac; /* Non-integer division factor */ 
    bool istat; /* Return value of CLI read function */

    sf_init(argc,argv);
    in = sf_input("in");
    inp_type = sf_gettype(in);
    n = sf_filesize(in);

    istat = sf_getint ("i_fac",&i_fac);
    
    if (inp_type == SF_FLOAT) {
        float_inp = true;
        nbuf = BUFSIZ/sizeof(float);
        fbuf = sf_floatalloc(nbuf);
        if (!istat) {
            if (!sf_getfloat("f_fac",&f_fac)) \
                sf_error("Need either i_fac= or f_fac=");
        }
        else f_fac = i_fac;
    }
    else if (inp_type == SF_INT) {
        int_inp = true;
        nbuf = BUFSIZ/sizeof(int);
        ibuf = sf_intalloc(nbuf);
        if (!istat) sf_error("Need i_fac=");
    }
    else sf_error("Need float or int input");

    /* Duplicating boilerplate code to avoid conditionals inside loops */

    if (float_inp) {

        for (; n > 0; n -= nbuf) {

            if (nbuf > n) nbuf = n;

            sf_floatread(fbuf, nbuf, in);

            for (i=0; i < nbuf; i++) {
                glob_i++;
                if (0 != fmod(fbuf[i],f_fac)) \
                    sf_error("Mismatch at element %d", glob_i);
            }
        }
    }

    else if (int_inp) {

        for (; n > 0; n -= nbuf) {

            if (nbuf > n) nbuf = n;

            sf_intread(ibuf, nbuf, in);

            for (i=0; i < nbuf; i++) {
                glob_i++;
                if (0 != (ibuf[i] % i_fac )) \
                    sf_error("Mismatch at element %d", glob_i);
            }
        }
    }


    exit(0);
}
