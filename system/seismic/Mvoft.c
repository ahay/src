/* V(t) function for a linear V(Z) profile.
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

int main(int argc, char** argv)
{
    int n1, n2, i1, i2;
    float* vint=NULL;
    float o1, d1, v0, alpha, t;
    sf_file in=NULL, out=NULL;

    sf_init (argc, argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint (in,"n1",&n1)) sf_error ("No n1= in input\n");
    if (!sf_histint (in,"n2",&n2)) sf_error ("No n2= in input\n");

    if (!sf_histfloat (in,"o1",&o1)) o1=0.;
    if (!sf_histfloat (in,"d1",&d1)) sf_error ("No d1= in input\n");

    if (!sf_getfloat ("v0",&v0)) v0 = 1.5;
    /* initial velocity */
    if (!sf_getfloat ("alpha",&alpha)) alpha = 0.5;
    /* velocity gradient */

    vint = sf_floatalloc(n1);
    for (i1=0; i1 < n1; i1++) {
	t = alpha*(o1+i1*d1);
	if (t > 0.) {
	    vint[i1] = v0 * sqrtf((expf(t)-1.)/t);
	} else {
	    vint[i1] = v0;
	}
    }
    for (i2=0; i2 < n2; i2++) {
	sf_floatwrite(vint,n1,out);
    }

    exit (0);
}
