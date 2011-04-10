/* 2-D synthetic model for multiple-arrival generation.

From Bill Symes.
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

#include "symes.h"

int main(int argc, char* argv[])
{
    int nx, nz, ix, iz;
    float* trace;
    float dx, dz, x[2];
    sf_file mod;

    sf_init (argc,argv);
    mod = sf_output("out");

    if (!sf_getint("nx",&nx)) nx=400; 
    /* horizontal dimension */
    if (!sf_getint("nz",&nz)) nz=800; 
    /* vertical dimension */

    dx = 1./(nx-1);
    dz = 2./(nz-1);

    sf_putint   (mod,"n1",nz); 
    sf_putfloat (mod,"d1",dz); 
    sf_putfloat (mod,"o1",0.);
    sf_putint   (mod,"n2",nx); 
    sf_putfloat (mod,"d2",dx); 
    sf_putfloat (mod,"o2",0.);
    sf_setformat (mod,"native_float");

    trace = sf_floatalloc(nz);

    for (ix = 0; ix < nx; ix++) {
	x[1] = ix*dx;
	for (iz = 0; iz < nz; iz++) {
	    x[0] = iz*dz;
	    trace[iz] = 1./sqrtf(symes_vel(NULL,x));
	}
	sf_floatwrite(trace,nz,mod);
    }

    exit (0);
}

/* 	$Id$	 */
