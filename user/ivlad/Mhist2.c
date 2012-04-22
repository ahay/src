/* Compute a 2-D histogram of integer- or float-valued input data
The output grid is not centered on the bins; it marks their "left edge".
I.e., the first sample holds the number of values between o1 and o1+d1*/
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

int main(int argc, char* argv[]) 
{
    int i;      /* Counter over input */
    int n;      /* Total number of values in input */
    int i1, i2, i12; /* Counters over output */
    int n1, n2; /* Output axes lengths */
    int n12;    /* = n1*n2 */
    int nbuf;   /* Number of elements read at one time */
    int   *hist=NULL; /* Output histogram */
    int   *ibuf1=NULL, *ibuf2=NULL; /* Input buffer for reading integers */
    float *fbuf1=NULL, *fbuf2=NULL; /* Input buffer for reading floats */
    float o1, d1, o2, d2; /* Output axes origin and sampling */
    sf_file in=NULL, inp2=NULL, out=NULL; /* Input, output files */
    sf_datatype inp_type; /* Input data type */

    sf_init(argc,argv);

    in   = sf_input("in");
    inp2 = sf_input("inp2");
    out  = sf_output("out");
    
    inp_type = sf_gettype(in);

    if(inp_type != sf_gettype(inp2)) sf_error("Input data types must match");

    if (inp_type == SF_FLOAT) {
        nbuf = BUFSIZ/sizeof(float);
        fbuf1 = sf_floatalloc(nbuf);
        fbuf2 = sf_floatalloc(nbuf);
    } else if (inp_type == SF_INT) {
        nbuf = BUFSIZ/sizeof(int);
        ibuf1 = sf_intalloc(nbuf);
        ibuf2 = sf_intalloc(nbuf);
    } else {
	nbuf = 0;
	sf_error("Need float or int input");
    }

    n = sf_filesize(in);

    if (!sf_getint("n1",&n1)) sf_error("Need n1=");
    /* number of histogram samples in dimension 1 */
    if (!sf_getfloat("o1",&o1)) sf_error("Need o1=");
    /* histogram origin for dimension 1 */
    if (!sf_getfloat("d1",&d1)) sf_error("Need d1=");
    /* histogram sampling for dimension 1 */
    if (!sf_getint("n2",&n2)) sf_error("Need n2=");
    /* number of histogram samples in dimension 2 */
    if (!sf_getfloat("o2",&o2)) sf_error("Need o2=");
    /* histogram origin for dimension 2 */
    if (!sf_getfloat("d2",&d2)) sf_error("Need d2=");
    /* histogram sampling for dimension 2 */

    sf_settype(out,SF_INT);
    sf_putint(out,"n1",n1);
    sf_putint(out,"n2",n2);
    sf_putfloat(out,"o1",o1);
    sf_putfloat(out,"o2",o2);
    sf_putfloat(out,"d1",d1);
    sf_putfloat(out,"d2",d2);

    /* If input n3... are defined, replace them with 1's in output */ 
       
    for (i=2; i < SF_MAX_DIM; i++) {
        int this_n;
        char key[5];
        sprintf(key,"n%d",i+1);
        if (!sf_histint(in,key,&this_n)) break;
 	    if (this_n > 1) sf_putint(out,key,1);
    }

    n12 = n1 * n2;

    hist = sf_intalloc(n12);
    for (i12=0; i12 < n12; i12++) {
        hist[i12]=0;
    }

    /* Duplicating boilerplate code to avoid conditionals inside loops */ 

    if (inp_type == SF_FLOAT) {

        for (; n > 0; n -= nbuf) {

            if (nbuf > n) nbuf = n;

            sf_floatread(fbuf1, nbuf, in);
            sf_floatread(fbuf2, nbuf, inp2);

            for (i=0; i < nbuf; i++) {

                i1 = (int) floorf((fbuf1[i]-o1)/d1);
                i2 = (int) floorf((fbuf2[i]-o2)/d2);
                if (i1 >= 0 && i1 < n1 && i2 >= 0 && i2 < n2) hist[i1+i2*n1]++;
            }
        }
    }
    
    else if (inp_type == SF_INT) {

        for (; n > 0; n -= nbuf) {

            if (nbuf > n) nbuf = n;                                        
                                                                           
            sf_intread(ibuf1, nbuf, in);                                    
            sf_intread(ibuf2, nbuf, inp2);
                                                                           
            for (i=0; i < nbuf; i++) {                                    
 
                i1 = (int) floorf((ibuf1[i]-o1)/d1);                        
                i2 = (int) floorf((ibuf2[i]-o2)/d2);                       
                if (i1 >= 0 && i1 < n1 && i2 >= 0 && i2 < n2) hist[i1+i2*n1]++;
            }                                                              
        }
    }

    sf_intwrite(hist, n12, out);
    
    exit(0);
}
