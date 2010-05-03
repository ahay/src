/* 2D pad and interleave traces.

Each initial trace is followed by "jump" zero traces.
*/
/*
Copyright (C) 2010 University of Texas at Austin

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

int main(int argc, char* argv[])
{
    int n1,n2,n3, jump, i1,i2,i3, j;
    float d, *pp=NULL, *zero=NULL, *one=NULL;
    sf_file in=NULL, out=NULL, mask=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    mask = sf_output("mask");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_getint("jump",&jump)) jump=2;
    /* aliasing */

    if (n2 > 1) {
	sf_putint(out,"n2",n2*jump);
	sf_putint(mask,"n2",n2*jump);
	if (sf_histfloat(in,"d2",&d)) {
	    sf_putfloat(out,"d2",d/jump);
	    sf_putfloat(mask,"d2",d/jump);
	}
    }

    pp = sf_floatalloc(n1);
    zero = sf_floatalloc(n1);
    one = sf_floatalloc(n1);

    for (i1=0; i1 < n1; i1++) {
	zero[i1] = 0.;
	one[i1] = 1.;
    }

    for (i3=0; i3 < n3; i3++) {
	sf_warning("slice %d of %d;",i3+1,n3);
	for (i2=0; i2 < n2; i2++) {
	    sf_floatread(pp,n1,in);
	    sf_floatwrite(pp,n1,out);
	    sf_floatwrite(one,n1,mask);
	    if (n2 >1) {
		for (j=0; j < jump-1; j++) {
		    sf_floatwrite(zero,n1,out);
		    sf_floatwrite(zero,n1,mask);
		}
	    }
	}
    }
    sf_warning(".");
    exit(0);
}
/* 	$Id$	 */
