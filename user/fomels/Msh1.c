/* Generate 1D shifts for regularized regression. */
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
    int n1, i2, n2, s1, ns;
    float *inp, **out;
    sf_file in, shift;

    sf_init(argc,argv);
    in = sf_input("in");
    shift = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getint("s1",&s1)) sf_error("Need s1=");
    /* shifts in samples */
    
    ns = 2*s1+1;
    sf_putint(shift,"n2", ns);
    sf_shiftdim(in, shift, 2);

    inp = sf_floatalloc(n1);
    out = sf_floatalloc2(n1,ns);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(inp,n1,in);
	shifts1(s1, n1, inp, out);
	sf_floatwrite(out[0],n1*ns,shift);
    }

    exit(0);
}
