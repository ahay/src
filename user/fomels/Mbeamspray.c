/* 2-D beam spraying. */
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

int main(int argc, char* argv[])
{
    int n1, nc, nd, n3, i3, nb, n1d, id, ic, i1;
    float *dense, *a, *p, *c, d2;
    sf_file in, out, dip, cur;

    sf_init(argc,argv);
    in = sf_input("in");
    dip = sf_input("dip");
    cur = sf_input("cur");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nc)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_histfloat(in,"d2",&d2)) d2=1.;

    nd = (nc-1)*nb+1;
    sf_putint(out,"n2",nd);
    sf_putfloat(out,"d2",d2/nb);

    n1d = n1*nd;

    dense = sf_floatalloc(n1d);
    a = sf_floatalloc(n1);
    p = sf_floatalloc(n1);
    c = sf_floatalloc(n1);
 
    for (i3=0; i3 < n3; i3++) {
	for (id=0; id < n1d; id++) {
	    dense[id] = 0.;
	}
	for (ic=0; ic < nc; ic++) {
	    sf_floatread(a,n1,in);
	    sf_floatread(p,n1,dip);
	    sf_floatread(c,n1,cur);

	    for (i1=0; i1 < n1; i1++) {
	    }
	}
	sf_floatwrite(dense,n1d,out);
    }

    exit(0);
}

