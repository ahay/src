/* 3D Coherence cube */

/*
  Copyright (C) 2012 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <rsf.h>

#include "coh1.h"


int main(int argc, char* argv[])
{
	sf_file in, out;
	int n1, n2, n3, n4, lag2, lag3, nw;
	int i3, i4;
	float **u1;

    sf_init(argc, argv);

    in  = sf_input("in");	/* 3D data set */
    out = sf_output("out");	/* 3D coherence cube */

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");
	n4 = sf_leftsize(in, 3);

    if (!sf_getint("nw", &nw) ) nw=5;
    /* half window size for coherence */
    if (!sf_getint("lag2",&lag2) ) lag2=3;
    /* maximal time lag on 2nd axis */
    if (!sf_getint("lag3",&lag3) ) lag3=3;
    /* maximal time lag on 3rd axis */

    u1 = sf_floatalloc2(n1, n2);

	coh1_init(nw, n1, n2, lag2, lag3);

	for(i4=0; i4<n4; i4++)
	{
		sf_floatread(u1[0], n1*n2, in);
		coh1_init2(u1);
		for(i3=1; i3<n3; i3++)
		{
			sf_floatread(u1[0], n1*n2, in);
			coh1(u1);
			sf_floatwrite(u1[0], n1*n2, out);
		}
		sf_floatwrite(u1[0], n1*n2, out);
		sf_warning("%d of %d", i4, n4);
	}
	coh1_close();
    free(u1[0]);
	free(u1);
    return 0;
}

