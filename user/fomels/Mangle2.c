/* Another illustration of angle gathers.
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
    int ix, iy, nx, ny;
    float dx, dy, zx, zy, x, y, *trace;
    sf_file angle;
    
    sf_init(argc,argv);
    angle=sf_output("out");
    sf_setformat(angle,"native_float");

    if (!sf_getint("nx",&nx)) nx=451;
    if (!sf_getint("ny",&ny)) ny=451;

    if (!sf_getfloat("dx",&dx)) dx=0.1;
    if (!sf_getfloat("dy",&dy)) dy=0.1;

    if (!sf_getfloat("zx",&zx)) zx=0.;
    if (!sf_getfloat("zy",&zy)) zy=0.;

    zx = tanf(SF_PI*zx/180.);
    zy = tanf(SF_PI*zy/180.);

    sf_putint(angle,"n1",2*nx-1);
    sf_putfloat(angle,"o1",-(nx-1)*dx);
    sf_putfloat(angle,"d1",dx);
    sf_putstring(angle,"label1","In-line Offset Slope");
    sf_putstring(angle,"unit1","degrees");

    sf_putint(angle,"n2",2*ny-1);
    sf_putfloat(angle,"o2",-(ny-1)*dy);
    sf_putfloat(angle,"d2",dy);
    sf_putstring(angle,"label2","Cross-line Offset Slope");
    sf_putstring(angle,"unit2","degrees");

    trace = sf_floatalloc(2*nx-1);
    
    for (iy=-ny+1; iy < ny; iy++) {
	y = tanf(iy*dy*SF_PI/180.);
	for (ix=-nx+1; ix < nx; ix++) {
	    x = tanf(ix*dx*SF_PI/180.);
	    x = x*x*(1.+zx*zx) + 2.*x*y*zx*zy + y*y*(1.+zy*zy);
	    x /= (1.+zx*zx+zy*zy);
	    if (x > 0.) {
		trace[ix+nx-1] = atanf(sqrtf(x)) * 180./SF_PI;
	    } else {
		trace[ix+nx-1] = -1.;
	    }
	}
	sf_floatwrite(trace,2*nx-1,angle);
    }

    exit(0);
}

/* 	$Id: Mangle2.c 7107 2011-04-10 02:04:14Z ivlad $	 */
