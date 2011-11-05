/* Compute 2-D Riesz transform. */
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

#include <rsf.h>

#include "riesz.h"

int main (int argc, char* argv[])
{
    int nx,ny,n2, i2, n;
    float **data, **datax, **datay, c;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&ny)) sf_error("No n2= in input");
    n2 = sf_leftsize(in,2);

    sf_shiftdim(in, out, 3);
    sf_putint(out, "n3", 2);

    data = sf_floatalloc2(nx,ny);
    datax = sf_floatalloc2(nx,ny);
    datay = sf_floatalloc2(nx,ny);

    if (!sf_getint("order",&n)) n=10;
    /* Hilbert transformer order */
    if (!sf_getfloat("ref",&c)) c=1.;
    /* Hilbert transformer reference (0.5 < ref <= 1) */

    riesz_init(n, c);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data[0],nx*ny,in);

	riesz(nx,ny,data,datax,datay);

	sf_floatwrite(datax[0],nx*ny,out);
	sf_floatwrite(datay[0],nx*ny,out);
    }

    exit(0);
}

