/* Analytical point-source traveltime in a linear slowness squared model. */
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

int main (int argc, char* argv[])
{
    int n, n1, n2, i1, i2;
    float d, g, s, v0, v, x, z, g2, v2;
    float **vel, **time[2];
    sf_file out;

    sf_init (argc,argv);
    out = sf_output("out");

    if (!sf_getint("n",&n)) sf_error("Need n=");
    /* number of samples */
    if (!sf_getfloat("g",&g)) g = -6.;    
    /* slowness squared gradient */
    g *= 0.5;
    g2 = g*g;
    if (!sf_getfloat("v0",&v0)) v0 = 4;
    /* initial slowness squared */
    if (!sf_getfloat("s",&s)) s = 0.5;
    /* shot location at the surface */
    
    d = 0.5/n;
    n1 = n+1;
    n2 = 3*n+1;

    sf_putint(out,"n1",n1); sf_putfloat(out,"d1",d); sf_putfloat(out,"o1",0.);
    sf_putint(out,"n2",n2); sf_putfloat(out,"d2",d); sf_putfloat(out,"o2",0.);
    sf_putint(out,"n3",3);
    sf_setformat(out,"native_float");

    vel = sf_floatalloc2 (n1,n2);
    time[0] = sf_floatalloc2 (n1,n2);
    time[1] = sf_floatalloc2 (n1,n2);

    for (i2 = 0; i2 < n2; i2++) {
	x = i2*d - s;
	x = x*x;
	for (i1 = 0; i1 < n1; i1++) {
	    z = i1*d;
	    v = v0 + g*z;
	    vel[i2][i1] = 1./sqrtf(v + g*z);
	    z = z*z+x;
	    if (v*v - g2*z >= 0.) {
		v2 = sqrt(v*v-g2*z);
		time[0][i2][i1] = sqrtf(2.*z/(v+v2))*(2.*v+v2)/3.;
		time[1][i2][i1] = sqrtf(2.*(v+v2)/g2)*(2.*v-v2)/3.;
	    } else {
		time[0][i2][i1] = -1.;
		time[1][i2][i1] = -1.;
	    }
	}
    }

    sf_floatwrite(vel[0], n1*n2,out);
    sf_floatwrite(time[0][0],n1*n2,out);
    sf_floatwrite(time[1][0],n1*n2,out);

    exit (0);
}

/* 	$Id$	 */
