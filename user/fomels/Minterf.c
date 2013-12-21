/* Create an interferometric matrix */
/*
  Copyright (C) 2013 University of Texas at Austin
   
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

#include "fft1.h"

int main(int argc, char* argv[])
{
    int n1, n2, nt, nw, i1, i2, i3;
    float *trace;
    sf_complex **fft, *ftrace;
    sf_file inp, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1 in input");
    if (!sf_histint(inp,"n2",&n2)) n2=1;
    sf_putint(out,"n3",n2);

    nt = fft1_init(n1,&nw);
    trace = sf_floatalloc(nt);
    fft = sf_complexalloc2(nw,n2);
    ftrace = sf_complexalloc(nw);

    ifft1_allocate(ftrace);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(trace,n1,inp);
	for (i1=n1; i1 < nt; i1++) {
	    trace[i1] = 0.0f;
	}
	fft1(trace,fft[i2]);
    }

    for (i3=0; i3 < n2; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < nw; i1++) {
		ftrace[i1] = fft[i3][i1]*conjf(fft[i2][i1]);
	    }
	    ifft1(trace,ftrace);
	    sf_floatwrite(trace,n1,out);
	}
    }

    exit(0);
}
