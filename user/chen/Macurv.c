/* Azimuth CURVature */

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
#include "acurv.h"
#include "runtime.h"

int main(int argc, char*argv[])
{
	sf_file in, out;
	int n1, n2, n3, rect[3], order, n12, nazim;
	int i1, i3;
	float a, **u1, *v;
	char *interp;

	sf_init(argc, argv);

	in  = sf_input("in");
	out = sf_output("out");
	
	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

	if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
	if (!sf_histint(in, "n2", &n2)) sf_error("No n2= in input");
	if (!sf_histint(in, "n3", &n3)) n3=1;

	if (!sf_getint("rect1",&rect[0])) rect[0]=1;
	/* smoothness on 1st axis */
	if (!sf_getint("rect2",&rect[1])) rect[1]=1;
	/* smoothness on 2nd axis */
	if (!sf_getint("rect3",&rect[2])) rect[2]=1;
	/* smoothness on 3rd axis */
	if (!sf_getint("order",&order)) order=2;
	/* approximating order of finite difference */
	if (!sf_getint("nazmuth",&nazim)) nazim=10;
	/* azimuth number */
	if ((interp=sf_getstring("interp"))==NULL) interp="maxflat";
	/* interpolation method: maxflat lagrange bspline */

	n12 =n1*n2;
	sf_shiftdim(in,out,3);
	sf_putint(out, "n3", nazim);
	sf_putint(out, "o3", 0);
	sf_putint(out, "d3", 1);
	v = sf_floatalloc(n12);
	u1 = sf_floatalloc2(n12, nazim);

	runtime_init(nazim*n12*sizeof(float));
	acurv_init(n1, n2, n3, order, interp, rect, nazim);
//  time sequences for recursive operators:
//          0   od   od+rc        n3  n3+od  n3+od+rc
//  i3      |----|-----|-----------|----|-----|
//  read    |----------------------|
//  acurv        |----------------------|
//  den          |----------------------|
//  write              |----------------------|

	for(i3=0; i3<n3+order+rect[2]; i3++)
	{
		if(i3<n3)
			sf_floatread(v, n12, in);

		for(i1=0; i1<n12; i1++) u1[0][i1] = v[i1];
		acurv(u1);

		if(i3>=order+rect[2])
		{
			sf_floatwrite(u1[0], n12*nazim, out);
			a = runtime(1);
			sf_warning("write %d of %d, %f MB/sec;", i3-order-rect[2], n3, a);
		}
		else sf_warning("read %d of %d;", i3, order+rect[2]);
	}

	acurv_close();
	free(u1[0]);
	free(u1);
	free(v);

	return 0;
}



