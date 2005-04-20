/* 2-D smoothing by triangle directional shaping. */
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

#include "trisl.h"

int main(int argc, char* argv[])
{
    bool adj;
    int n1, n2, n12, n3, i3, rect1, rect2;
    float *input, *smooth, **slope;
    sf_file in, out, dip;

    sf_init(argc,argv);
    in = sf_input("in");
    dip = sf_input("dip");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    n12 = n1*n2;

    if (!sf_getint("rect1",&rect1)) rect1=3;
    if (!sf_getint("rect2",&rect2)) rect2=3;
    /* smoothing radius */

    if (!sf_getbool("adj",&adj)) adj=false;

    trisl_init(n1, n2, rect1, rect2);

    input = sf_floatalloc(n12);
    smooth = sf_floatalloc(n12);
    slope = sf_floatalloc2(n1,n2);
 
    for (i3=0; i3 < n3; i3++) {
	sf_floatread(input,n12,in);
	sf_floatread(slope[0],n12,dip);

	trisl_set(slope);

	if (adj) {
	    trisl_lop(true,false,n12,n12,smooth,input);
	} else {
	    trisl_lop(false,false,n12,n12,input,smooth);
	}

	sf_floatwrite(smooth,n12,out);
    }

    exit(0);
}

/* 	$Id$	 */
