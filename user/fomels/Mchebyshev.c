/* Testing Chebyshev interpolation */
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

#include "chebyshev.h"

int main (int argc, char* argv[])
{
    int nc, i1, n1;
    float o1, d1, x;
    float *data, *intp;
    sf_file inp, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&nc)) sf_error("No n1= in input");

    if (!sf_getint("n1",&n1)) n1=1; /* number of output points */
    if (!sf_getfloat("o1",&o1)) o1=0.0f; /* output origin */
    if (!sf_getfloat("d1",&d1)) d1=1.0f; /* output sampling */

    sf_putint(out,"n1",n1);
    sf_putfloat(out,"o1",o1);
    sf_putfloat(out,"d1",d1);

    data = sf_floatalloc(nc);
    intp = sf_floatalloc(n1);

    chebyshev_init(nc,o1,o1+d1*(n1-1));

    sf_floatread(data,nc,inp);
    chebyshev_set(data);

    for (i1=0; i1 < n1; i1++) {
	x = o1 + i1*d1;
	intp[i1] =  chebyshev(x);
    }

    sf_floatwrite(intp,n1,out);

    exit(0);
}
