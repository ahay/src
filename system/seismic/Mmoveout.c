/* Put spikes at an arbitrary moveout */
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
#include "stretch4.h"

int main(int argc, char* argv[])
{
    map4 mo;
    int n1, i2, n2, iw, nw, ns;
    float o1, d1, t, eps;
    float *trace, *move, *str, *amp;
    sf_file out, warp;

    sf_init(argc,argv);
    warp = sf_input("in");
    out = sf_output("out");

    if (!sf_getint("n1",&n1)) sf_error("Need n1="); /* time samples */
    if (!sf_getfloat("d1",&d1)) d1=1.; /* time sampling */
    if (!sf_getfloat("o1",&o1)) o1=0.; /* time origin */

    sf_shiftdim(warp, out, 1);
	
    sf_putint(out,"n1",n1);
    sf_putfloat(out,"d1",d1);
    sf_putfloat(out,"o1",o1);
    sf_putstring(out,"label1","Time");
    sf_putstring(out,"unit1","s");

    n2 = sf_filesize(warp);

    if (!sf_getfloat("eps",&eps)) eps=0.1;
    /* stretch regularization */

    if (!sf_getint("nw",&nw)) nw=10;
    /* wavelet length */
    ns = 2*nw+1;

    move = sf_floatalloc(n2);
    sf_floatread(move,n2,warp);

    trace = sf_floatalloc(n1);
    str = sf_floatalloc(ns);
    amp = sf_floatalloc(ns);

    mo = stretch4_init (n1, o1, d1, ns, eps);

    for (i2=0; i2 < n2; i2++) {
	t = move[i2];
	str[nw] = t;
	amp[nw] = 1.0;

	for (iw=0; iw < nw; iw++) {
	    str[iw] = t - (nw-iw)*d1;
	    amp[iw] = 0.0;

	    str[nw+iw+1] = t+iw*d1;
	    amp[nw+iw+1] = 0.0;
	}

	stretch4_define (mo,str);
	stretch4_apply (mo,amp,trace);
	sf_floatwrite (trace,n1,out);
    }

    exit(0);
}
