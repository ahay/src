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
	int n1, n2, n3, n4, lag1, lag2, nw;
	int i3, i4;
	float **u1;
	bool twod, verb;

	sf_init(argc, argv);

	in  = sf_input("in");	/* 3D data set */
	out = sf_output("out");	/* 3D coherence cube */

	if (!sf_getint("nw", &nw) ) nw=5;
	/* half window size for coherence */
	if (!sf_getint("lag1",&lag1) ) lag1=3;
	/* maximal time lag on 2nd axis */
	if (!sf_getint("lag2",&lag2) ) lag2=3;
	/* maximal time lag on 3rd axis */
	if (!sf_getbool("twod",&twod) ) twod=true;
	/* y: only twod coherence */
	if (!sf_getbool("verb",&verb) ) verb=true;
	/* verbosity */

	if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
	if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
	if (!sf_histint(in,"n3",&n3)) n3=1; 
	if(n3==1)  twod=true;

	u1 = sf_floatalloc2(n1, n2);

	coh1_init(nw, n1, n2, lag1, lag2);

	if(twod) n4 = sf_leftsize(in, 2);
	else 	n4 = sf_leftsize(in, 3);

	for(i4=0; i4<n4; i4++)
	{
		sf_floatread(u1[0], n1*n2, in);
		if(twod)
		{
			coh1_2d(u1);
			if(verb) sf_warning("%d of %d;", i4, n4);
		}else{
			coh1_3d_init2d(u1);
			for(i3=1; i3<n3; i3++)
			{
				sf_floatread(u1[0], n1*n2, in);
				coh1_3d(u1);
				sf_floatwrite(u1[0], n1*n2, out);
				if(verb) sf_warning("%d of %d;", i3+i4*n3, n3*n4);
			}
		}
		sf_floatwrite(u1[0], n1*n2, out);
	}
	coh1_close();
	free(u1[0]);
	free(u1);
	return 0;
}

