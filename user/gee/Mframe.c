/* Create a frame for binning.
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

static float frame(float x, float y, float base, int nd, float** xyz)
{
    int id;
    float wt, sum, m;

    m = 0.;
    sum = 0.;
    for (id=0; id < nd; id++) {
	wt = hypotf(x-xyz[id][0],y-xyz[id][1]);
	wt = 1/(wt*wt);
	wt *= wt; /* inverse distance to the 4th power */
	m += (xyz[id][2]-base)*wt;
	sum += wt;
    }
    m /= sum;

    return m;
}

int main(int argc, char* argv[]) {
    int ix, iy, nx, ny, nd, three;
    float **xyz, *mm, dx, dy, x0, y0, x, y, base;
    sf_file in, xyzs, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    xyzs = sf_input("xyz");

    if (!sf_histint(in,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&ny)) sf_error("No n2= in input");

    if (!sf_histfloat(in,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&dy)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o1",&x0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"o2",&y0)) sf_error("No o2= in input");

    mm = sf_floatalloc(nx);

    if (!sf_histint(xyzs,"n1",&three) || three != 3) 
	sf_error("Need n1=3 in xyz");
    if (!sf_histint(xyzs,"n2",&nd)) sf_error("No n2= in xyz");

    xyz = sf_floatalloc2(3,nd);

    sf_floatread(xyz[0],3*nd,xyzs);

    if (!sf_getfloat("base",&base)) base=0.;
    /* base to be subtracted from z */

    for (iy=0; iy < ny; iy++) {
	y = y0+iy*dy;

	sf_floatread(mm,nx,in);

	if (0==iy || ny-1==iy) {
	    for (ix=0; ix < nx; ix++) {
		x = x0+ix*dx;
		mm[ix] = frame(x,y,base,nd,xyz);
	    }
	} else {
	    mm[0]    = frame(x0,          y,base,nd,xyz);
	    mm[nx-1] = frame(x0+(nx-1)*dx,y,base,nd,xyz);
	}

	sf_floatwrite(mm,nx,out);
    }

    exit(0);
}

/* 	$Id$	 */
