/* orientation estimation by structural gradient tensor */

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
#include "rgradient.h"
#include "vecfilt.h"

int main(int argc, char*argv[])
{
	sf_file in, wgt, out;
	int n1, n2, n3, rect[3], order;
	int i1, i3;
	float *u1, *w, *v1, *v2;
	char *interp;

	sf_init(argc, argv);

	in  = sf_input("in");
	wgt  = sf_input("weight");
	out = sf_output("out");

	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

	if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
	if (!sf_histint(in, "n2", &n2)) sf_error("No n2= in input");
	if (!sf_histint(in, "n3", &n3)) n3=1;

	if (!sf_getint("rect1",&rect[0])) rect[0]=0;
	/* smoothness on 1st axis */
	if (!sf_getint("rect2",&rect[1])) rect[1]=0;
	/* smoothness on 2nd axis */
	if (!sf_getint("rect3",&rect[2])) rect[2]=0;
	/* smoothness on 3rd axis */
	if (!sf_getint("order",&order)) order=4;
	/* approximating order of finite difference */
	if ((interp=sf_getstring("interp"))==NULL) interp="maxflat";
	/* interpolation method: maxflat lagrange bspline */

	rgradient_init(interp, order, n1, n2);
	vecfilt_init(3, n1, n2, rect);

	u1 = sf_floatalloc(n1*n2);
	v1 = sf_floatalloc(n1*n2*3);
	v2 = sf_floatalloc(n1*n2*3);

	if(wgt) w  = sf_floatalloc(n1*n2);
	else w = NULL;

	for(i3=0; i3<n3; i3++)
	{	
		sf_floatread(u1, n1*n2, in);
		if(wgt) sf_floatread(w, n1*n2, wgt);
		rgradient(u1, v1);
		vecfilt(v1, v2, w);
		for(i1=0; i1<n1*n2; i1++)
		u1[i1] = atanf(fabs(v2[i1*3]) / 
			sqrt(v2[i1*3+1]*v2[i1*3+1]+v2[i1*3+2]+v2[i1*3+2]));
		sf_floatwrite(u1, n1*n2, out);
	}

	rgradient_close();
	vecfilt_close();

	free(u1);
	free(v1);
	free(v2);

	if(w) free(w);

	return 0;
}



