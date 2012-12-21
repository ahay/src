/* 2d filter bank  */

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
#include "fbank.h"


int main(int argc, char*argv[])
{
	sf_file in, out;
	int nf, n1, n2, n3, m, n;
	int i3;
	float **wav, ****fb;
	char *interp;

	sf_init(argc, argv);

	in  = sf_input("in");
	out = sf_output("out");

	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

	if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
	if (!sf_histint(in, "n2", &n2)) sf_error("No n2= in input");
	n3 = sf_leftsize(in, 2);

	sf_shiftdim2(in,out,2);

	if(!sf_getint("m", &m)) m=1;
	/* b[-m, ... ,n] */
	if(!sf_getint("n", &n)) n=1;
	/* b[-m, ... ,n] */
	if ((interp=sf_getstring("interp"))==NULL) interp="maxflat";
	/* interpolation method: maxflat lagrange bspline */

	nf = m+n+1;

	wav = sf_floatalloc2(n1,n2);
	fb  = sf_floatalloc4(n1, n2, nf, nf);

	sf_putint(out, "n3", nf);
	sf_putint(out, "n4", nf);
	sf_putfloat(out, "o3", 0);
	sf_putfloat(out, "d3", 1);
	sf_putfloat(out, "o4", 0);
	sf_putfloat(out, "d4", 1);


	fbank_init(m, n, interp);


	for(i3=0; i3<n3; i3++)
	{
		sf_floatread(wav[0], n1*n2, in);
		fbank2(n1, n2, wav, fb);
		sf_floatwrite(fb[0][0][0], n1*n2*nf*nf, out);
	}

	fbank_close();
	return 0;
}



