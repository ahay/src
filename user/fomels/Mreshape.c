/* Non-stationary spectral balancing. */
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

#include "reshape.h"

int main(int argc, char* argv[])
{
    int n1, n2, n3, i2, i3, iw, n, nfft, dim;
    float o1, d1, *data, *data2=NULL, m1, a1, m2, a2;
    char key[5];
    sf_file in, in2, ma, ma2, out, out2=NULL;

    sf_init(argc,argv);
    in = sf_input("in");

    if (NULL != sf_getstring("in2")) {
	/* optional second input file */
	in2 = sf_input("in2");
    } else {
	in2 = NULL;
    }
    
    ma = sf_input("ma");
    ma2 = sf_input("ma2");

    out = sf_output("out");

    if (NULL != in2) out2 = sf_output("out2");

    if ((SF_FLOAT != sf_gettype (in)) ||
	(NULL != in2 && 
	 SF_FLOAT != sf_gettype (in2))) sf_error("Need float data");

    if (!sf_histint(in,"n1",&n1)) n1=1;

    if (NULL != in2 && sf_histint(in2,"n1",&n) && n != n1)
	sf_error("Size mismatch in in2: %d != %d",n,n1);

    if (!sf_getint("dim",&dim)) dim=1;
    /* data dimensionality */
    sprintf(key,"n%d",dim+1);

    if (!sf_histint(in,key,&n3)) n3=1;

    if (NULL != in2 && sf_histint(in2,key,&n) && n != n3)
	sf_error("Size mismatch in in2 [%s]: %d != %d",key,n,n3);

    n2 = sf_leftsize(in,1);
    n2 /= n3;

    /* n3 is the number of windows, n2xn1 is the window size */

    if (!sf_histfloat(in,"d1",&d1)) d1=0.004;
    if (!sf_histfloat(in,"o1",&o1)) o1=0.;
    
    /* determine frequency sampling (for real to complex FFT) */
    nfft = 2*kiss_fft_next_fast_size((n1+1)/2);
    reshape_init (nfft,d1);

    if (!sf_histint(ma,"n1",&n) || n != 2)
	sf_error("Wrong n1= in ma");
    if (!sf_histint(ma2,"n1",&n) || n != 2)
	sf_error("Wrong n1= in ma2");
    if (!sf_histint(ma,"n2",&n) || n != n3)
	sf_error("Wrong n2= in ma");
    if (!sf_histint(ma2,"n2",&n) || n != n3)
	sf_error("Wrong n2= in ma2");
    
    data = sf_floatalloc(nfft);
    if (NULL != in2) data2 = sf_floatalloc(nfft);

    for (i3=0; i3 < n3; i3++) { /* loop over windows */
	sf_floatread(&m1,1,ma);
	sf_floatread(&a1,1,ma);
	sf_floatread(&m2,1,ma2);
	sf_floatread(&a2,1,ma2);
	
	for (i2=0; i2 < n2; i2++) { /* loop over traces in a window */
	    sf_floatread(data,n1,in);
	    if (NULL != in2) sf_floatread(data2,n1,in2);
	
	    /* Fourier transform */
	    if (m1 > m2) {
		for (iw=n1; iw < nfft; iw++) {
		    data[iw] = 0.;
		}
		reshape (m1, a1, m2, a2, data);
	    } else if (NULL != in2) {
		for (iw=n1; iw < nfft; iw++) {
		    data2[iw] = 0.;
		}
		reshape (m2, a2, m1, a1, data2);
	    }
	    
	    sf_floatwrite(data,n1,out);
	    if (NULL != in2) sf_floatwrite(data2,n1,out2);
	}
    }

    exit(0);
}

/* 	$Id$	 */
