/* Convolution with a trapezoidal filter. */
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

#include "trapez.h"

int main(int argc, char* argv[])
{
    int i, n1, n2, i2;
    float d1, freq[4], *trace=NULL;
    sf_file in=NULL, out=NULL;

    sf_init(argc,argv);
    in  = sf_input ( "in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getfloats("frequency",freq,4)) {
	/* frequencies (in Hz), default: (0.1,0.15,0.45,0.5)*Nyquist */
	freq[0] = 0.1*0.5;
	freq[1] = 0.15*0.5;
	freq[2] = 0.45*0.5;
	freq[3] = 0.5*0.5;
    } else {
	if (!sf_histfloat(in,"d1",&d1)) d1=1.;
	for (i=0; i < 4; i++) {
	    freq[i] *= d1;
	}
    }

    trace = sf_floatalloc(n1);
    trapez_init(2*kiss_fft_next_fast_size((n1+1)/2),freq);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(trace,n1,in);
	sf_freqfilt(n1,trace);
	sf_floatwrite(trace,n1,out);
    }

    exit(0);
}
