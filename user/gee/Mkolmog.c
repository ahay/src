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

int main(int argc, char* argv[]) 
{
    int i1, n1, nfft, nw, lag, shift;
    float *trace;
    bool spec;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("spec",&spec)) spec=false;
    /* if y, the input is spectrum squared; n, time-domain signal */

    if (!sf_getint("lag",&lag)) lag=0;
    /* lag for asymmetric part */

    if (!sf_getint("shift",&shift)) shift=0;
    /* time shift */

    if (spec) { 
	if (!sf_histint(in,"n1",&nw)) sf_error("No n1= in input");
	nfft = 2*(nw-1);

	sf_putint(out,"n1",nfft);
	sf_settype(out,SF_FLOAT);

	trace = sf_floatalloc(nfft);

	sf_floatread(trace,nw,in);
	kolmog_init(nfft,lag,shift);
	kolmog2(trace);
	
	sf_floatwrite(trace,nfft,out);
    } else { /* input signal */
	if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");

	/* determine wavenumber sampling (for real to complex FFT) */
	nfft = 2*kiss_fft_next_fast_size((n1+1)/2);
	trace = sf_floatalloc(nfft);
	
	sf_floatread(trace,n1,in);

	for (i1=n1; i1 < nfft; i1++) {
	    trace[i1]=0.;
	}

	kolmog_init(nfft,lag,shift);
	kolmog(trace);

	sf_floatwrite(trace,n1,out);
    }

    exit(0);
}

/* 	$Id: Mkolmog.c 10702 2013-07-25 21:37:48Z sfomel $	 */
