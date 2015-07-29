/* Remove random shot gathers from a 2-D dataset. */
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

int main(int argc, char* argv[])
{
    int n2, n3, i2, i3, is, **known=NULL;
    float *chance=NULL, perc;
    sf_file in=NULL, mask=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    mask = sf_output("out");

    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");

    sf_putint(mask,"n1",n2);
    sf_putint(mask,"n2",n3);
    sf_putint(mask,"n3",1);
    sf_settype(mask,SF_INT);

    if (!sf_getfloat("perc",&perc)) perc=0.75;
    /* how many shots to remove */

    known = sf_intalloc2(n2,n3);
    chance = sf_floatalloc(n2+n3);

    init_genrand(2003);
    sf_random (n2+n3,chance);

    for (i3=0; i3 < n3; i3++) { /* half-offset */
	for (i2=0; i2 < n2; i2++) { /* midpoint */
	    is = i2 - i3 + n3-1; /* shot */
	    known[i3][i2] = (chance[is] > perc);
	}
    }
    sf_intwrite (known[0],n2*n3,mask);

    exit(0);
}
