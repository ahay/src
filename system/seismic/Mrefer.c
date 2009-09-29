/* Subtract a reference from a grid. */
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
#include "fint1.h"

int main(int argc, char* argv[])
{
    int n1, n2, n3, m1, i1, i2, i3, i;
    float *trace=NULL, **datum=NULL, o1, d1, f, shift;
    fint1 fnt;
    sf_file in=NULL, out=NULL, ref=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    m1 = 2*n1-1;
    trace = sf_floatalloc(m1);
    datum = sf_floatalloc2(n2,n3);

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_histfloat(in,"o1",&o1)) o1=0.;

    sf_putint(out,"n1",m1);
    sf_putfloat(out,"o1",-(n1-1)*d1);

    ref = sf_input("ref");
    sf_floatread(datum[0],n2*n3,ref);

    fnt = fint1_init (5, n1, 0);

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    sf_floatread(trace,n1,in);
	    fint1_set(fnt,trace);
	    shift = datum[i3][i2] - o1;
	    for (i1=0; i1 < m1; i1++) {
		f = ((i1-n1+1)*d1+shift)/d1;
		i = f;
		if (i >=0 && i < n1-1) {
		    trace[i1] = fint1_apply (fnt, i, f-i, false);
		} else {
		    trace[i1] = 0.;
		}
	    }
	    sf_floatwrite(trace,m1,out);
	}
    }

    exit(0);
}
