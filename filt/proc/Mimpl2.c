/* 2-D anisotropic diffusion.
*/
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

#include "impl2.h"

int main(int argc, char* argv[])
{
    int n1, n2, n3, i3;
    float r1, r2, tau, pclip, **dat;
    bool up;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_getfloat("rect1",&r1)) sf_error("Need rect1=");
    /* vertical smoothing */
    if (!sf_getfloat("rect2",&r2)) sf_error("Need rect2=");
    /* horizontal smoothing */
    if (!sf_getfloat("tau",&tau)) tau=0.1;
    /* smoothing time */
    if (!sf_getfloat("pclip",&pclip)) pclip=50.;
    /* percentage clip for the gradient */
    if (!sf_getbool ("up",&up)) up=false;
    /* smoothing style */

    impl2_init (r1, r2, n1, n2, tau, pclip, up);

    dat = sf_floatalloc2(n1,n2);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread (dat[0],n1*n2,in);
	impl2_apply (dat,true,false);
	sf_floatwrite (dat[0],n1*n2,out);
    }

    sf_close();
    exit(0);
}

/* 	$Id: Mimpl2.c,v 1.5 2004/06/25 18:08:42 fomels Exp $	 */

