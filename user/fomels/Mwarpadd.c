/* Add a perturbation to the warping function.
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

#include <math.h>

#include <rsf.h> 

int main(int argc, char* argv[])
{ 
    int i1, n1, i2, n2, m2, i3, n3, order, i;
    float *first, *second, o1, d1, x, f, f1, t;
    sf_eno ent;
    sf_file in, sum, add;

    sf_init (argc, argv);
    in = sf_input("in");
    add = sf_input("add");
    sum = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if(!sf_histint(in,"n1",&n1)) sf_error ("No n1= in input");
    if(!sf_histfloat(in,"d1",&d1)) sf_error ("No d1= in input");
    if(!sf_histfloat(in,"o1",&o1)) o1 = 0.;

    if(!sf_histint(in,"n2",&m2)) m2 = 1;
    n3 = sf_leftsize(in,2);

    if(!sf_getint("accuracy",&order)) {
	/* Interpolation accuracy order */
	order = 2;
    } else if (order < 1 || order > 4) {
	sf_error ("accuracy must be between 1 and 4");
    }
    order *= 2;

    if (!sf_getint("m1",&n2)) n2=n1*2;
    /* Trace pading */

    ent = sf_eno_init(order,n2);

    first = sf_floatalloc (n1); 
    second = sf_floatalloc (n2);

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < m2; i2++) {
	    sf_floatread(first,n1,add);
	    sf_floatread(second,n1,in);
	    
	    for (i1=n1; i1 < n2; i1++) {
		second[i1] = (7.*second[i1-1]-5.*second[i1-2]+second[i1-3])/3.;
	    }

	    sf_eno_set (ent,second);
	    
	    for (i1=0; i1 < n1; i1++) {
		t = (o1+i1*d1);
		x = (first[i1] + t - o1)/d1;
		i = floorf(x); x -= i;
		sf_eno_apply (ent, i, x, &f, &f1, FUNC);
		first[i1] += f;
	    }
	    
	    sf_floatwrite(first,n1,sum);
	}
    }

    exit (0);
}

