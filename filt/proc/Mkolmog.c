/* Kolmogoroff spectral factorization. */
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

#include "kolmog.h"

int main(int argc, char* argv[]) {
    int i1, n1, nfft, nw;
    float *trace;
    float complex *fft;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX == sf_gettype(in)) { /* input complex spectrum */
	if (!sf_histint(in,"n1",&nw)) sf_error("No n1= in input");
	nfft = 2*(nw-1);

	sf_putint(out,"n1",nfft);
	sf_settype(out,SF_FLOAT);

	trace = sf_floatalloc(nfft);
	fft = sf_complexalloc(nw);

	sf_complexread(fft,nw,in);
	kolmog2(nfft,nw,trace,fft);
	sf_floatwrite(trace,nfft,out);
    } else { /* input signal */
	if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");

	/* determine wavenumber sampling (for real to complex FFT) */
	nfft = sf_npfar(n1);
	trace = sf_floatalloc(nfft);
	
	sf_floatread(trace,n1,in);
	for (i1=n1; i1 < nfft; i1++) {
	    trace[i1] = 0.;
	}
	kolmog(nfft, trace);
	sf_floatwrite(trace,n1,out);
    }
    
    exit(0);
}

/* 	$Id: Mkolmog.c,v 1.3 2004/07/02 11:54:47 fomels Exp $	 */
