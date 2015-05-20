/* Find minimum in 2-D */
/*
  Copyright (C) 2015 University of Texas at Austin
  
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

#include "findmin2.h"

int main(int argc, char* argv[])
{
    int i1, i2, i3, n1, n2, n3, gate1, gate2, ib1, ie1, ib2, ie2, l1, l2, m1, m2, k1, k2;
    float **data, **prob, c, d, f, xy[2], fmin[3];
    sf_file inp, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(inp,2);

    data = sf_floatalloc2(n1,n2);

    sf_putint(out,"n1",3);
    sf_putint(out,"n2",1);

    if (!sf_getint("gate1",&gate1)) gate1=3;
    if (!sf_getint("gate2",&gate2)) gate2=3;
    /* picking gate */

    prob = sf_floatalloc2(gate1,gate2);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(data[0],n1*n2,inp);

	fmin[0] = SF_HUGE;

	for (i2=0; i2 < n2; i2++) {
	    ib2 = SF_MAX(i2-gate2,-1);
	    ie2 = SF_MIN(i2+gate2,n2);

	    for (i1=0; i1 < n1; i1++) {
		ib1 = SF_MAX(i1-gate1,-1);
		ie1 = SF_MIN(i1+gate1,n1);
			
		c = SF_HUGE;
		l1 = -1;
		l2 = -1;

		for (k2=ib2+1; k2 < ie2; k2++) {
		    m2 = k2-ib2-1;
		    for (k1=ib1+1; k1 < ie1; k1++) {
			m1 = k1-ib1-1;
			
			d = data[k2][k1];
			
			if (d < c) {
			    c = d;
			    l1 = m1;
			    l2 = m2;
			}

			prob[m2][m1]=d;
		    }
		}

		/* local minimum */
		f = find_minimum(l1,ie1-ib1-1,ib1+1,
				 l2,ie2-ib2-1,ib2+1,
				 c,prob,xy);
		
		if (f < fmin[0]) {
		    fmin[0] = f;
		    fmin[1] = xy[0];
		    fmin[2] = xy[1];
		}
	    }
	}

	sf_floatwrite(fmin,3,out);
    }
    
    exit(0);					
}
