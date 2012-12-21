/* fault detection  */

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

int main(int argc, char*argv[])
{
	int n1, n2, i1, i2;
	sf_file in, out;
	float *u1, *v, *p, *q, lamda, sigma;

	sf_init(argc, argv);

	in = sf_input("in");
	out = sf_output("out");

	if(!sf_getfloat("lamda", &lamda)) lamda=0.9;
	if(!sf_getfloat("sigma", &sigma)) sigma=1.0;
	if(!sf_histint(in, "n1", &n1)) sf_error("n1 needed in input");
	n2 = sf_leftsize(in, 1);

	u1 = sf_floatalloc(n1);
	v = sf_floatalloc(n1);
	p = sf_floatalloc(n1);
	q = sf_floatalloc(n1);

	for(i2=0; i2<n2; i2++)
	{
		sf_floatread(u1, n1, in);
		v[0] = u1[0];
		for(i1=1; i1<n1; i1++)
		{
			p[i1] = (u1[i1] - v[i1-1])/sigma;
			p[i1] *= p[i1];
			p[i1] = expf(-p[i1]);
			v[i1] = lamda * p[i1];
			v[i1] = (1-v[i1])*u1[i1] + v[i1]*v[i1-1];
		}
		v[n1-1] = u1[n1-1];
		for(i1=n1-2; i1>=0; i1--)
		{
			q[i1] = (u1[i1] - v[i1+1])/sigma;
			q[i1] *= q[i1];
			q[i1] = expf(-q[i1]);
			v[i1] = lamda * q[i1];
			v[i1] = (1-v[i1])*u1[i1] + v[i1]*v[i1+1];
		}
		for(i1=1; i1<n1-1; i1++) v[i1] = (p[i1] * q[i1]);
		v[0] = v[1]; v[n1-1] = v[n1-2];
		sf_floatwrite(v, n1, out);
	}

	free(u1);
	free(p);
	free(q);
	free(v);
}

