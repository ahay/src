/* Shot interpolation. */
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
    sf_complex *s=NULL, **ss=NULL; /* a[3]; */
    sf_file in=NULL, out=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");

    if (!sf_histint(in,"n1",&nh)) sf_error("No n1= in input");
    /* offsets */
    if (!sf_histint(in,"n2",&ns)) sf_error("No n2= in input");
    /* shots */

    if (!sf_histint(in,"n3",&nw)) sf_error("No n3= in input");
    /* frequency */

    if (!sf_histfloat(in,"o1",&h0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&dh)) sf_error("No d1= in input");

    if (!sf_histfloat(in,"d2",&ds)) sf_error("No d2= in input");

    if (!sf_histfloat(in,"o3",&w0)) sf_error("No o3= in input");
    if (!sf_histfloat(in,"d3",&dw)) sf_error("No d3= in input");

    if (!sf_getfloat("eps",&eps)) eps=0.1;
    /* regularization parameter */
    eps *= eps;

    if (!sf_getbool("positive",&sign)) sign=true;
    /* initial offset orientation */

    sf_putint(out,"n2",2*ns-1);
    sf_putfloat(out,"d2",0.5*ds);
    /* make the shot spacing denser */

    dw *= 2.*SF_PI;
    w0 *= 2.*SF_PI;
    /* convert Hertz to radian */

    ss = sf_complexalloc2(nh,ns);
    s = sf_complexalloc(nh);
    /* allocate space for shots at one frequency slice */

    shotfill_init(nh,h0,dh,sign? 0.5*ds: -0.5*ds, eps);

    for (iw=0; iw < nw; iw++) {
	/* loop over frequency slices */
	w = w0 + iw*dw;

	sf_complexread (ss[0],nh*ns,in);

	if(fabsf(w) < dw) { /* dc */
	    /* write out zeroes */
	    for (ih=0; ih < nh; ih++) {
		s[ih] = sf_cmplx(0.,0.);
	    }
	    sf_complexwrite(s,nh,out);
	    for (is=1; is < ns; is++) {
		sf_complexwrite(s,nh,out);
		sf_complexwrite(s,nh,out);
	    }
	    continue;
	}

	shotfill_define(w);
	/* set coefficients */

	sf_complexwrite(ss[0],nh,out);
	/* write first shot */

	for (is=1; is < ns; is++) {
	    /* loop over shots */

	    shotfill_apply(ss[is-1],ss[is],s);
	    /* insert shot between is-1 and is */

	    sf_complexwrite(s,nh,out);
	    sf_complexwrite(ss[is],nh,out);
	} /* s */
    } /* w */

    exit(0);
}
