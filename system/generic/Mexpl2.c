/* 2-D anisotropic diffusion by box cascade. */
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

#include "expl2.h"

int main(int argc, char* argv[])
{
    int n1, n2, n12, n3, i3, rect, ic, nc;
    float pclip, **dat;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    n12 = n1*n2;

    if (!sf_getint("cycle",&nc)) nc=1;
    /* number of cycles */
    if (!sf_getint("rect",&rect)) sf_error("Need rect=");
    /* vertical smoothing */
    if (!sf_getfloat("pclip",&pclip)) pclip=50.;
    /* percentage clip for the gradient */

    expl2_init (rect, n1, n2, pclip);

    dat = sf_floatalloc2(n1,n2);

    for (i3=0; i3 < n3; i3++) {
	sf_warning("slice %d of %d;",i3,n3);

	sf_floatread (dat[0],n12,in);
	for (ic=0; ic < nc; ic++) {
	    expl2_apply (dat);
	}
	sf_floatwrite (dat[0],n12,out);
    }
    sf_warning(".");

    exit(0);
}
