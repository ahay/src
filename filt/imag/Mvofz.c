/* Analytical traveltime in a linear V(z) model. */
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
    float g, s, v0, v, x, z, a, d1, d2, o1, o2;
    float **vel, **time;
    bool intime;
    sf_file out;

    sf_init (argc,argv);
    out = sf_output("out");

    if (!sf_getint("n1",&n1)) {
	/* vertical samples */
	if (!sf_getint("n",&n)) sf_error("Need n=");
	/* number of samples */
	n1 = n+1;
    }
    if (!sf_getint("n2",&n2)) {
	/* horizontal samples */
	if (!sf_getint("n",&n)) sf_error("Need n=");
	/* number of samples */
	n2 = 2*n+1;
    }
    if (!sf_getfloat("g",&g)) g = 1.;
    /* velocity gradient */
    if (!sf_getfloat("v0",&v0)) v0 = 0.5;
    /* initial velocity */
    if (!sf_getfloat("s",&s)) s = 0.5;
    /* shot location at the surface */
    
    if (!sf_getfloat("d1",&d1)) d1 = 0.5/(n1-1);
    /* vertical sampling */
    if (!sf_getfloat("d2",&d2)) d2 = 1./(n2-1);
    /* horizontal sampling */
    if (!sf_getfloat("o1",&o1)) o1 = 0.;
    /* vertical origin */
    if (!sf_getfloat("o2",&o2)) o2 = 0.;
    /* horizontal origin */

    if (!sf_getbool("intime",&intime)) intime = false;
    /* if in vertical time coordinates */

    a = 0.5*g*g/v0;

    sf_putint(out,"n1",n1); sf_putfloat(out,"d1",d1); sf_putfloat(out,"o1",0.);
    sf_putint(out,"n2",n2); sf_putfloat(out,"d2",d2); sf_putfloat(out,"o2",0.);
    sf_putint(out,"n3",2);
    sf_setformat(out,"native_float");

    vel = sf_floatalloc2 (n1,n2);
    time = sf_floatalloc2 (n1,n2);

    for (i2 = 0; i2 < n2; i2++) {
	x = o2+i2*d2 - s;
	x = x*x;
	for (i1 = 0; i1 < n1; i1++) {
	    z = o1+i1*d1;
	    if (intime) z = v0*(expf(g*z)-1.)/g;

	    v = v0 + g*z;
	    z = 1. + a*(z*z+x)/v;
	    vel[i2][i1] = v;
	    time[i2][i1] = fabsf(logf(z + sqrtf(z*z-1.))/g);
	}
    }

    sf_floatwrite(vel[0], n1*n2,out);
    sf_floatwrite(time[0],n1*n2,out);
    
    exit (0);
}

/* 	$Id$	 */
