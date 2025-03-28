/* Generate 2D shifts for regularized regression. */
/*
  Copyright (C) 2025 University of Texas at Austin
  
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

#include "shifts.h"

int main(int argc, char* argv[])
{
    int n1, n2, i3, n3, s1, s2, ns;
    float **inp, ***out;
    sf_file in, shift;

    sf_init(argc,argv);
    in = sf_input("in");
    shift = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_getint("s1",&s1)) sf_error("Need s1=");
    if (!sf_getint("s2",&s2)) sf_error("Need s2=");
    /* shifts in samples */
    
    ns = (2*s1+1)*(s2+1);
    sf_putint(shift,"n3", ns);
    sf_shiftdim(in, shift, 3);

    inp = sf_floatalloc2(n1,n2);
    out = sf_floatalloc3(n1,n2,ns);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(inp[0],n1*n2,in);
	shifts2(s1, s2, n1, n2, inp, out);
	sf_floatwrite(out[0][0],n1*n2*ns,shift);
    }

    exit(0);
}
