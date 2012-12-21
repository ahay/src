/* Gabor transform by linear phase filter bank */

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
#include "dgauss.h"

int main(int argc,char*argv[])
{
	sf_file in, out;
	int i1, i2, j2, i3, n1, n2, n3, nf;
	float **u1, **u2, sigma, d2;

	sf_init(argc,argv);

	in  = sf_input("in");
	out = sf_output("out");

	n3 = sf_leftsize(in, 2);

	if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input file");
	if (!sf_histint(in, "n2", &n2)) sf_error("No n2= in input file");

	if (!sf_getfloat("sigma", &sigma)) sigma=1.0;
	/* sigma */
	if (!sf_getint("nf", &nf)) nf=100;
	/* frequency samples [0, 0.5] */

	d2 = SF_PI/(nf-1);

	sf_putint(out, "n2", nf);
	sf_putfloat(out, "d2", 0.5/(nf-1));
	dgauss_init(n1, 1.0/sigma);

	u1 = sf_floatalloc2(n1, n2);
	u2 = sf_floatalloc2(n1, nf);

	for (i3=0; i3<n3; i3++)
	{
		sf_floatread(u1[0], n1*n2, in);
		for(i1=0; i1<n1; i1++)
		{
			for (j2=0; j2<nf; j2++)
			for (i2=0; i2<n2; i2++)
			{
				u2[j2][i1] = u1[i2][i1] * dgauss_val(i2, d2*j2) 
					* powf(0.5/SF_PI, i2);
			}
		}
		sf_floatwrite(u2[0], n1*nf, out);
	}

	free(u1[0]);
	free(u2[0]);
	free(u1);
	free(u2);
	dgauss_close();
	exit(0);
}

