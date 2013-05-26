/* horizon extraction  */

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
#include "sinterp.h"

int main(int argc, char*argv[])
{
	int n1, n2, n3,n4, i2, i3, i4, m2, m3;
	sf_file in, hor, out;
	float **u, **v, **h, o1, d1;
	char *interp;
	sinterp intp;

	sf_init(argc, argv);

	in = sf_input("in");
	hor = sf_input("horizon");
	out = sf_output("out");

	if(!sf_histint(in, "n1", &n1)) sf_error("n1 needed in input");
	if(!sf_histfloat(in, "o1", &o1)) o1=0.0;
	if(!sf_histfloat(in, "d1", &d1)) d1=1.0;
	if(!sf_histint(in, "n2", &n2)) sf_error("n2 needed in input");
	if(!sf_histint(in, "n3", &n3)) n3=1; 
	if(!sf_histint(in, "n4", &n4)) n4=1; 
	if(!sf_histint(hor, "n1", &m2)) sf_error("n1 needed in horizon file");
	if(!sf_histint(hor, "n2", &m3)) m3=1; 

	if((interp=sf_getstring("interp")) ==NULL ) interp="linear";
	/*< interpolation method: nearest, linear >*/

	intp = sinterp_c2f(interp);
	if(n2!=m2 || n3!= m3) sf_error("horizon file not match");

	sf_unshiftdim(in,out,1);
	u = sf_floatalloc2(n1, n2);
	v = sf_floatalloc2(n2, n3);
	h = sf_floatalloc2(n2, n3);

	sf_floatread(h[0], n3*n2, hor);
	for(i4=0; i4<n4; i4++)
	{
		for(i3=0; i3<n3; i3++)
		{
			sf_floatread(u[0], n1*n2, in);
			for(i2=0; i2<n2; i2++)
			v[i3][i2] = intp(u[i2], (h[i3][i2]-o1)/d1, n1);
		}
		sf_floatwrite(*v, n2*n3, out);
	}

	free(*u);
	free(u);
	free(*v);
	free(v);
	free(*h);
	free(h);
}

