/* Compute a histogram of integer- or float-valued input data.

The output grid is not centered on the bins; it marks their "left edge".
I.e., the first sample holds the number of values between o1 and o1+d1. 

February 2015 program of the month:
http://ahay.org/blog/2015/03/01/program-of-the-month-sfhistogram/
*/
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
    int i;  /* Counter over input */
    int n;  /* Total number of values in input */
    int i1; /* Counter over output */
    int n1; /* Output axis length */
    int nbuf; /* Number of elements read at one time */
    int   *hist; /* Output histogram */
    int   *ibuf=NULL; /* Input buffer for reading integers */
    float *fbuf=NULL; /* Input buffer for reading floats */
    float o1, d1;     /* Output axis origin and sampling */
    sf_file in, out; /* Input, output files */
    sf_datatype inp_type; /* Input data type */

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    inp_type = sf_gettype(in);

    if (inp_type == SF_FLOAT) {
        nbuf = BUFSIZ/sizeof(float);
        fbuf = sf_floatalloc(nbuf);
    } else if (inp_type == SF_INT) {
        nbuf = BUFSIZ/sizeof(int);
        ibuf = sf_intalloc(nbuf);
    } else {
	nbuf = 0;
	sf_error("Need float or int input");
    }

    n = sf_filesize(in);

    if (!sf_getint("n1",&n1)) sf_error("Need n1=");
    /* number of histogram samples */
    if (!sf_getfloat("o1",&o1)) sf_error("Need o1=");
    /* histogram origin */
    if (!sf_getfloat("d1",&d1)) sf_error("Need d1=");
    /* histogram sampling */

    sf_settype(out,SF_INT);
    sf_putint(out,"n1",n1);
    sf_putfloat(out,"o1",o1);
    sf_putfloat(out,"d1",d1);

    /* If input n2, n3... are defined, replace them with 1's in output */

    for (i=1; i < SF_MAX_DIM; i++) {
        int this_n;
        char key[5];
        sprintf(key,"n%d",i+1);
        if (!sf_histint(in,key,&this_n)) break;
 	    if (this_n > 1) sf_putint(out,key,1);
    }

    hist = sf_intalloc(n1);
    for (i1=0; i1 < n1; i1++) {
        hist[i1]=0;
    }

    /* Duplicating boilerplate code to avoid conditionals inside loops */

    switch (inp_type) {
	case SF_FLOAT:
	    for (; n > 0; n -= nbuf) {
		
		if (nbuf > n) nbuf = n;
		
		sf_floatread(fbuf, nbuf, in);
		
		for (i=0; i < nbuf; i++) {
		    i1 = (int) floorf((fbuf[i]-o1)/d1);
		    if (i1 >= 0 && i1 < n1) hist[i1]++;
		}
	    }
	    break;
	case SF_INT:
	    for (; n > 0; n -= nbuf) {
		
		if (nbuf > n) nbuf = n;                                        
		
		sf_intread(ibuf, nbuf, in);                                    
		
		for (i=0; i < nbuf; i++) {                                     
		    i1 = (int) floorf((ibuf[i]-o1)/d1);                        
		    if (i1 >= 0 && i1 < n1) hist[i1]++;                        
		}                                                              
	    }
	    break;
	default:
	    sf_error("Need float or int input");
	    break;
    }

    sf_intwrite(hist,n1,out);

    exit(0);
}
