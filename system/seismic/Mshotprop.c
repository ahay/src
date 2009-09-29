/* Shot propagation. */
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

#include <math.h>
#include <rsf.h>
#include "shotfill.h"

int main(int argc, char* argv[])
{
    int ns, nh, nw, iw, ih, is;
    bool sign;
    float ds, h0, dh, w0, dw, w, eps;
    sf_complex *s=NULL, *s2=NULL;
    sf_file in=NULL, out=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");

    if (!sf_histint(in,"n1",&nh)) sf_error("No n1= in input");
    /* offsets */
    if (!sf_histint(in,"n3",&nw)) sf_error("No n2= in input");
    /* frequency */

    if (!sf_histfloat(in,"o1",&h0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&dh)) sf_error("No d1= in input");

    if (!sf_histfloat(in,"o3",&w0)) sf_error("No o3= in input");
    if (!sf_histfloat(in,"d3",&dw)) sf_error("No d3= in input");

    if (!sf_getint("ns",&ns)) sf_error("Need ns=");
    /* number of shots */
    if (!sf_getfloat("ds",&ds)) sf_error("Need ds=");
    /* shot sampling */

    sf_putint(out,"n2",ns);
    sf_putfloat(out,"d2",ds);

    if (!sf_getfloat("eps",&eps)) eps=0.1;
    /* regularization parameter */
    eps *= eps;

    if (!sf_getbool("positive",&sign)) sign=true;
    /* initial offset orientation */

    dw *= 2.*SF_PI;
    w0 *= 2.*SF_PI;
    /* convert Hertz to radian */

    s2 = sf_complexalloc(nh);
    s = sf_complexalloc(nh);
    /* allocate space for shots at one frequency slice */

    shotfill_init(nh,h0,dh,sign? ds: -ds, eps);

    for (iw=0; iw < nw; iw++) {
	/* loop over frequency slices */
	w = w0 + iw*dw;

	sf_complexread (s,nh,in);

	if(fabsf(w) < dw) { /* dc */
	    /* write out zeroes */
	    for (ih=0; ih < nh; ih++) {
		s2[ih] = sf_cmplx(0.,0.);
	    }
	    for (is=0; is < ns; is++) {
		sf_complexwrite(s2,nh,out);
	    }
	    continue;
	}

	shotprop_define(w);
	/* set coefficients */
	sf_complexwrite(s,nh,out);

	for (is=1; is < ns; is++) {
	    /* loop over shots */

	    shotprop_apply(s,s2);
	    /* extrapolate shot */
	    for (ih=0; ih < nh; ih++) {
		s[ih] = s2[ih];
	    }

	    sf_complexwrite(s,nh,out);
	} /* s */
    } /* w */

    exit(0);
}
