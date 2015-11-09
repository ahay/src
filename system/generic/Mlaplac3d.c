/* 3-D finite-difference Laplacian operation. */
/*
  Copyright (C) 2004 University of Texas at Austin
  Copyright (C) 2014 Colorado School of Mines

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

#include "laplac3.h"

int main(int argc, char* argv[])
{
    int n1,n2,n3, n123, n4,i4;
    float *pp,*qq;
    bool adj;
    sf_file in,out;

    sf_init(argc,argv);

    in  = sf_input(  "in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");
    n4 = sf_leftsize(in,3);
    n123 = n1*n2*n3;

    if (!sf_getbool("adj",&adj)) adj=false; /* adjoint flag */

    pp = sf_floatalloc(n123);
    qq = sf_floatalloc(n123);

    laplac3_init(n1,n2,n3);

    for (i4=0; i4<n4; i4++) {
	sf_floatread(pp,n123,in);
	if (adj) laplac3_lop ( true,false,n123,n123,qq,pp);
	else     laplac3_lop (false,false,n123,n123,pp,qq);
	sf_floatwrite(qq,n123,out);
    }

    exit(0);
}
