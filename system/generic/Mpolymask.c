/* Mask a polygon. */
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

#include "pnpoly.h"

int main(int argc, char* argv[])
{
    int i1, n1, i2, n2, nv, two, *trace;
    float o1, o2, d1, d2, x, y, **vert;
    sf_file inp, out, poly;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");

    if (!sf_histfloat(inp,"d1",&d1)) d1=1.;
    if (!sf_histfloat(inp,"d2",&d2)) d2=1.;
    
    if (!sf_histfloat(inp,"o1",&o1)) o1=0.;
    if (!sf_histfloat(inp,"o2",&o2)) o2=0.;

    poly = sf_input("poly"); 
    /* list of polygon vertices */

    if (SF_FLOAT != sf_gettype(poly)) sf_error("Need float type in poly");
    if (!sf_histint(poly,"n1",&two) || 2 != two)
	sf_error("Need n1=2 in poly");
    if (!sf_histint(poly,"n2",&nv))
	sf_error("No n2= in poly");
    
    trace = sf_intalloc(n1);
    vert = sf_floatalloc2(2,nv);

    sf_floatread(vert[0],2*nv,poly);

    sf_settype(out,SF_INT);

    for (i2=0; i2 < n2; i2++) {
	y = o2+i2*d2;
	for (i1=0; i1 < n1; i1++) {
	    x = o1+i1*d1;
	    trace[i1] = pnpoly(nv, vert, x, y);
	}
	sf_intwrite(trace,n1,out);
    }

    exit(0);
}
