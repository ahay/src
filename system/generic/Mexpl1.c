/* 1-D anisotropic diffusion by explicit cascade. */
/*
Copyright (C) 2019 University of Texas at Austin

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

#include "expl1.h"

int main(int argc, char* argv[])
{
    int n1, n2, i2, rect, ic, nc;
    float pclip, *dat=NULL;
    sf_file in=NULL, out=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getint("cycle",&nc)) nc=1;
    /* number of cycles */
    if (!sf_getint("rect1",&rect)) sf_error("Need rect1=");
    /* smoothing radius */
    if (!sf_getfloat("pclip",&pclip)) pclip=50.;
    /* percentage clip for the gradient */

    expl1_init (rect, n1, pclip);

    dat = sf_floatalloc(n1);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread (dat,n1,in);
	for (ic=0; ic < nc; ic++) {
	    expl1_apply (dat);
	}
	sf_floatwrite (dat,n1,out);
    }

    exit(0);
}
