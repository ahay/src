/* Analytical plane-wave traveltime in a linear slowness squared model. */
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
    float d, gz, gx, v0, v, x, z, z2, g2, v2;
    float **vel, **time;
    sf_file out;

    sf_init (argc,argv);
    out = sf_output("out");

    if (!sf_getint("n",&n)) sf_error("Need n=");
    /* number of samples */
    if (!sf_getfloat("gz",&gz)) gz = -6.;
    if (!sf_getfloat("gx",&gx)) gx = 2.;
    /* slowness squared gradient */
    gz *= 0.5;
    gx *= 0.5;
    g2 = 4*gx*gx+gz*gz;
    if (!sf_getfloat("v0",&v0)) v0 = 4;
    /* initial slowness squared */
    
    d = 0.5/n;
    n1 = n+1;
    n2 = 3*n+1;

    sf_putint(out,"n1",n1); sf_putfloat(out,"d1",d); sf_putfloat(out,"o1",0.);
    sf_putint(out,"n2",n2); sf_putfloat(out,"d2",d); sf_putfloat(out,"o2",0.);
    sf_putint(out,"n3",2);
    sf_setformat(out,"native_float");

    vel = sf_floatalloc2 (n1,n2);
    time = sf_floatalloc2 (n1,n2);

    for (i2 = 0; i2 < n2; i2++) {
	x = i2*d;
	for (i1 = 0; i1 < n1; i1++) {
	    z = i1*d;
	    z2 = g2*z*z;

	    v = v0 + gz*z+2*gx*x;
	    v2 = v*v;
	    vel[i2][i1] = 1./sqrtf(v + gz*z);
	    
	    
	    time[i2][i1] = z*(v2+z2/3.0)*
		sqrtf(2.0/(v*(v2+3.0*z2)+(v2-z2)*sqrtf(v2-z2)));
	}
    }

    sf_floatwrite(vel[0], n1*n2,out);
    sf_floatwrite(time[0],n1*n2,out);

    exit (0);
}

/* 	$Id$	 */
