/* Variogram from irregular 2-D data. */
/*
  Copyright (C) 2013 University of Texas at Austin

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

int main (int argc, char* argv[])
{
    int id, id2, nd, ix, nx, three, *fold;
    float *vari, **xy, distance, dx, diff;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    /* read data */

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float");
    if (!sf_histint(in,"n1",&three) || 3 != three) 
	sf_error("Need n1=3 in in");
    if (!sf_histint(in,"n2",&nd)) sf_error("Need n2=");

    xy = sf_floatalloc2(3,nd);
    sf_floatread(xy[0],3*nd,in);

    /* set output dimensions */

    if (!sf_getint ("nx",&nx)) sf_error("Need nx=");
    /* number of distance bins */
    if (!sf_getfloat("dx",&dx)) sf_error("Need dx=");
    /* distance sampling */

    sf_putint(out,"n1",nx);
    sf_putint(out,"n2",1);

    sf_putfloat (out,"o1",0.0f);
    sf_putfloat (out,"d1",dx);
 
    vari = sf_floatalloc(nx);
    fold = sf_intalloc(nx);

    for (ix=0; ix < nx; ix++) {
	fold[ix] = 0;
	vari[ix] = 0.0f;
    }

    /* loop through data pairs */
    for (id=0; id < nd; id++) {
	for (id2=id+1; id2 < nd; id2++) {
	    distance = hypotf(xy[id][0]-xy[id2][0],
			      xy[id][1]-xy[id2][1]);
	    ix = floorf(distance/dx);
	    if (ix < nx) {
		diff = xy[id][2]-xy[id2][2];
		vari[ix] += diff*diff;
		fold[ix]++;
	    }
	}
    }

    for (ix=0; ix < nx; ix++) {
	if (fold[ix] > 0) vari[ix] *= 0.5/fold[ix];
    }
    
    sf_floatwrite (vari,nx,out);
    exit(0);
}

