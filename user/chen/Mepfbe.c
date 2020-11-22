/* Bi-Exponential Edge Preserving Smoothing */

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
#include "beeps.h"


int main(int argc, char*argv[])
{
	sf_file in,  out;
	int i1, i2, i3, n1, n2, n3;
	float lamda, **u1, **u2=NULL, **u3=NULL, sigma;
	bool twod;
	kernel op;
	void *h1, *h2=NULL;
	char *str;

	sf_init(argc, argv);

	in = sf_input("in");
	out = sf_output("out");

	if(!sf_getbool("twod", &twod)) twod=false;
	/* y, 2D smoothing */
	if ((str=sf_getstring("kernel"))==NULL) str="gaussian";
	/* similarity: gaussian */
	if (!sf_getfloat("sigma", &sigma)) sigma=1.0;
	/* normalizing parameter */
	if(!sf_getfloat("lamda", &lamda)) lamda=0.8;
	/* lamda */

	op = kernel_c2f(str);

	if(!sf_histint(in, "n1", &n1)) sf_error("n1 needed in input");
	if(!sf_histint(in, "n2", &n2)) {n2=1; twod=false;}
	n3 = sf_leftsize(in, 2);

	u1 = sf_floatalloc2(n1, n2);
	h1 = beeps_init(n1, op, sigma, lamda);
	if(twod)
	{
		u2 = sf_floatalloc2(n1, n2);
		u3 = sf_floatalloc2(n1, n2);
		h2 = beeps_init(n2, op, sigma, lamda);
	}

	for(i3=0; i3<n3; i3++)
	{
		sf_floatread(u1[0], n1*n2, in);
		if(twod)
		{
			memcpy(u2[0], u1[0], n1*n2*sizeof(float));
			for(i2=0; i2<n2; i2++)	beeps(h1, u2[i2], 1);
			for(i1=0; i1<n2; i1++)	beeps(h2, u2[0]+i1, n1);
			memcpy(u3[0], u1[0], n1*n2*sizeof(float));
			for(i1=0; i1<n2; i1++)	beeps(h2, u3[0]+i1, n1);
			for(i2=0; i2<n2; i2++)	beeps(h1, u3[i2], 1);
			for(i1=0; i1<n1*n2; i1++)
				u1[0][i1] = 0.5*(u2[0][i1]+u3[0][i1]);
		} else
		for(i2=0; i2<n2; i2++)	beeps(h1, u1[i2], 1);
		sf_floatwrite(u1[0], n1*n2, out);
	}

	free(u1[0]);
	free(u1);
	if(twod)
	{
		free(u2[0]);
		free(u2);
		free(u3[0]);
		free(u3);
	}
	return 0;

}


