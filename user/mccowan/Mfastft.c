/* Fast Fourier Transform. */
/*
  Copyright (C) 2008 University of Texas at Austin

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

#include "fastfft.h"

int main (int argc, char *argv[])
{
    bool inv;
    int n1, i1, i2, n2, nt;
    float *real, *imag;
    sf_complex *data;
    sf_file in=NULL, out=NULL;

    sf_init(argc, argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, perform inverse transform */

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");

    if (!sf_histint  (in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    nt = fft_init(n1);
    data = sf_complexalloc(n1);
    real = sf_floatalloc(nt);
    imag = sf_floatalloc(nt);

    for (i2=0; i2 < n2; i2++) {
	sf_complexread(data,n1,in);

	for (i1=0; i1 < n1; i1++) {
	    real[i1] = crealf(data[i1]);
	    imag[i1] = cimagf(data[i1]);
	}

	for (i1=n1; i1 < nt; i1++) {
	    real[i1] = 0.;
	    imag[i1] = 0.;
	}

	fft(real,imag,inv);

	for (i1=0; i1 < n1; i1++) {
	    data[i1] = sf_cmplx(real[i1],imag[i1]);
	}

	sf_complexwrite(data,n1,out);
    }


    exit(0);
}
