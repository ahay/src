/* Creates a simple spherical surface. */
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
    int i1, i2, n1, n2;
    float o1, o2, d1, d2, x, y, *z;
    sf_file sphere;

    sf_init (argc,argv);
    sphere = sf_output("out");
    sf_setformat(sphere,"native_float");

    if (!sf_getint("n1",&n1)) n1=101; sf_putint(sphere,"n1",n1);
    if (!sf_getint("n2",&n2)) n2=101; sf_putint(sphere,"n2",n2);
    if (!sf_getfloat("o1",&o1)) o1=0.; sf_putfloat(sphere,"o1",o1);
    if (!sf_getfloat("o2",&o2)) o2=0.; sf_putfloat(sphere,"o2",o2);
    if (!sf_getfloat("d1",&d1)) d1=0.01; sf_putfloat(sphere,"d1",d1);
    if (!sf_getfloat("d2",&d2)) d2=0.01; sf_putfloat(sphere,"d2",d2);

    z = sf_floatalloc(n1);

    for (i2=0; i2 < n2; i2++) {
	y = o2+i2*d2;
	y = 0.25 - (y-0.5)*(y-0.5);
	for (i1 =0; i1 < n1; i1++) {
	    x = o1 + i1*d1;
	    x = y - (x-0.5)*(x-0.5);
	    z[i1] = (x > 0.)? sqrtf(x): 0.;
	}
	sf_floatwrite(z,n1,sphere);
    }

    exit(0);
}

