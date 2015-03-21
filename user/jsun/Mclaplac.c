/* 2-D finite-difference Laplacian operation for complex numbers. */
/*
  Copyright (C) 2015 University of Texas at Austin

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

#include "laplac2.h"

int main(int argc, char* argv[])
{
    int n1,n2,n12, n3,i3;
    sf_complex *pp,*qq;
    sf_file in,out;

    sf_init(argc,argv);

    in  = sf_input ( "in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    n12 = n1*n2;

    pp = sf_complexalloc(n12);
    qq = sf_complexalloc(n12);

    laplac2_init(n1,n2);

    for (i3=0; i3<n3; i3++) {
	sf_complexread(pp,n12,in);
        laplac2 (n12,n12,pp,qq);
	sf_complexwrite(qq,n12,out);
    }

    exit(0);
}
