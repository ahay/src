/* Bandpass filtering using erf function. */
/*
  Copyright (C) 2006 University of Texas at Austin
   
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
#include <rsf.h>

int main(int argc, char* argv[]) {
    int n1, nt, iw, nw, i2, n2;
    float fhi, flo, a, fi, d1, dw, f, *trace, *filt;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    n2 = sf_leftsize(in,1);
    
    nt = 2*kiss_fft_next_fast_size((n1+1)/2);
    nw = nt/2+1;
    sf_freqfilt_init(nt,nw);

    trace = sf_floatalloc(n1);
    filt = sf_floatalloc(nw);
    dw = 1./(nt*d1);

    if (!sf_getfloat("flo",&flo)) flo=-1.;
    /* low frequency in band */
    if (!sf_getfloat("fhi",&fhi)) fhi=-1.;
    /* high frequency in band */
    if (!sf_getfloat("a",&a)) a=0.25;
    /* filter sharpness */

    for (iw=0; iw < nw; iw++) {
	f = iw*dw;
	fi = 1.;
	if (flo > 0.) fi *= (1.-0.5*(erf(a*(f+flo))-erf(a*(f-flo))));
	if (fhi > 0.) fi *= 0.5*(erf(a*(f+fhi))-erf(a*(f-fhi)));
	filt[iw]=fi/nt;
    }

    sf_freqfilt_set(filt);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(trace,n1,in);
	sf_freqfilt(n1,trace);
	sf_floatwrite(trace,n1,out);
    }

    exit(0);
}
