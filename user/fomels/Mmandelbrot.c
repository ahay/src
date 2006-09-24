/* Generate Mandelbrot set. */
/*
  Copyright (C) 2006 University of Texas at Austin
   
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

int main(int argc, char* argv[]) 
{
    int i1, i2, n1, n2, *set, iter, niter;
    float xmax, ymax, x0, y0, dx, dy;
    double x, y, xc, yc, x2, y2;
    sf_file out;

    sf_init(argc,argv);
    out = sf_output("out");
    sf_setformat(out,"native_int");

    if (!sf_getint("n1",&n1)) n1=512;
    if (!sf_getint("n2",&n2)) n2=512;
    /* dimensions */

    sf_putint(out,"n1",n1);
    sf_putint(out,"n2",n2);

    if (!sf_getfloat("x0",&x0)) x0=-2.;
    if (!sf_getfloat("y0",&y0)) y0=-1.;
    /* set origin */

    sf_putfloat(out,"o1",x0);
    sf_putfloat(out,"o2",y0);

    if (!sf_getfloat("xmax",&xmax)) xmax=0.5;
    if (!sf_getfloat("ymax",&ymax)) ymax=1.;
    /* set maximum */

    if (!sf_getint("niter",&niter)) niter=1000;
    /* number of iterations */

    if (!sf_getfloat("dx",&dx)) dx = (xmax-x0)/(n1-1);
    if (!sf_getfloat("dy",&dy)) dy = (ymax-y0)/(n2-1);

    sf_putfloat(out,"d1",dx);
    sf_putfloat(out,"d2",dy);

    set = sf_intalloc(n1);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    x = xc = x0 + i1*dx;
	    y = yc = y0 + i2*dy;

	    x2 = x*x;
	    y2 = y*y;

	    for (iter=0; iter < niter; iter++) {
		if (x2 + y2 >= 4) break;

		y = 2*x*y + yc;
		x = x2 - y2 + xc;

		x2 = x*x;     
		y2 = y*y;
	    }

	    if (iter==niter) iter=0;
	    set[i1] = iter;
	}

	sf_intwrite(set,n1,out);
    }

    exit(0);
}
