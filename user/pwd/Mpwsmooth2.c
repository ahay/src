/* 2-D structure-enhancing filtering: two slopes. */
/*
  Copyright (C) 2014 University of Texas at Austin
  
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

#include "pwsmooth2.h"

int main(int argc, char* argv[])
{
    bool adj, verb;
    int n1, n2, n12, n3, i3, order, ns;
    float *input, *smooth, **slope1, **slope2, eps;
    sf_file in, out, dip;

    sf_init(argc,argv);
    in = sf_input("in");
    dip = sf_input("dip");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    n12 = n1*n2;

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (!sf_getint("ns",&ns)) ns=0;
    /* smoothing radius */

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */	

    input = sf_floatalloc(n12);
    smooth = sf_floatalloc(n12);
    slope1 = sf_floatalloc2(n1,n2);
    slope2 = sf_floatalloc2(n1,n2);

    pwsmooth2_init(ns, n1, n2, order, eps, slope1, slope2);

    for (i3=0; i3 < n3; i3++) {
	if (verb) sf_warning("slice %d of %d;",i3+1,n3);

	sf_floatread(input,n12,in);
	sf_floatread(slope1[0],n12,dip);
	sf_floatread(slope2[0],n12,dip);

	if (adj) {
	    pwsmooth2_lop(true,false,n12,n12,smooth,input);
	} else {
	    pwsmooth2_lop(false,false,n12,n12,input,smooth);
	}

	sf_floatwrite(smooth,n12,out);
    }
    if (verb) sf_warning(".");

    exit(0);
}


