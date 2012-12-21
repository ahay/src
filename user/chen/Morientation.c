/* orientation estimation by structural gradient tensor */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
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
#include "runtime.h"

int main(int argc, char*argv[])
{
	sf_file in, wgt, out;
	int n1, n2, n3, rect[3], order;
	int i1, i3;
	float *u1, *w, *v1, *v2, a;
	char *interp;

	sf_init(argc, argv);

	in  = sf_input("in");
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
	if(sf_getstring("weight")!=NULL) wgt  = sf_input("weight");
	else wgt = NULL;

	rgradient_init(interp, order, n1, n2);
	vecfilt_init(3, n1, n2, rect);

	u1 = sf_floatalloc(n1*n2);
	v1 = sf_floatalloc(n1*n2*3);
	v2 = sf_floatalloc(n1*n2*3);

	if(wgt) w  = sf_floatalloc(n1*n2);
	else w = NULL;

	runtime_init(n1*n2*sizeof(float));

//	time sequences for recursive operators:
//			0	od	 od+rc		  n3  n3+od	 n3+od+rc
//	i3		|----|-----|-----------|----|-----|
//	read	|----------------------|
//	u1=0						   |----|
//	grad	|----|----------------------|
//	v1=0								|-----|
//	weight		 |----------------------|
//	vecf		 |-----|----------------------|
//	write			   |----------------------|

	for(i3=0; i3<n3+order+rect[2]; i3++)
	{
		if(i3<n3) sf_floatread(u1, n1*n2, in);
		if(i3>=n3 && i3<n3+order) memset(u1, 0, n1*n2*sizeof(float));
		if(i3<n3+order)	rgradient(u1, v1);
		else memset(v1, 0, 3*n1*n2*sizeof(float));
		if(i3<order) continue;

		if(wgt && i3<n3+order) sf_floatread(w, n1*n2, wgt);
		vecfilt(v1, v2, w);
		if(i3<order+rect[2]) continue;

		for(i1=0; i1<n1*n2; i1++)
		{
			a = sqrt(v2[i1*3+1]*v2[i1*3+1]+v2[i1*3+2]*v2[i1*3+2]);
			u1[i1] = atan2(fabs(v2[i1*3]), a);
		}
		sf_floatwrite(u1, n1*n2, out);
		a = runtime(1);
		sf_warning("%d of %d, %f MB/sec;", i3, n3, a);
	}

	rgradient_close();
	vecfilt_close();

	free(u1);
	free(v1);
	free(v2);

	if(wgt) free(w);

	return 0;
}



