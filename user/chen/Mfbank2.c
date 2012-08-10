/* 2d filter bank  */

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
#include "fbank.h"


int main(int argc, char*argv[])
{
	sf_file in, out;
	int nf, n1, n2, n3, interp, nf2;
	int i3;
	float **wav, ****fb;

	sf_init(argc, argv);

	in  = sf_input("in");
	out = sf_output("out");

	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

	if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
	if (!sf_histint(in, "n2", &n2)) sf_error("No n2= in input");
	n3 = sf_leftsize(in, 2);

	sf_shiftdim2(in,out,2);

	if (!sf_getint("order",&nf)) nf=1;
	/* PWD filter order */
	
	if (!sf_getint("interp",&interp)) interp=0;
	/* interpolation method: 
	0: maxflat
	1: Lagrange 
	2: B-Spline */

	nf2 = 2*nf+1;
	wav = sf_floatalloc2(n1,n2);
	fb  = sf_floatalloc4(n1, n2, nf2, nf2);

	sf_putint(out, "n3", nf2);
	sf_putint(out, "n4", nf2);
	sf_putfloat(out, "o3", 0);
	sf_putfloat(out, "d3", 1);
	sf_putfloat(out, "o4", 0);
	sf_putfloat(out, "d4", 1);


	fbank_init(nf, interp);


	for(i3=0; i3<n3; i3++)
	{
		sf_floatread(wav[0], n1*n2, in);
		fbank2(n1, n2, wav, fb);
		sf_floatwrite(fb[0][0][0], n1*n2*nf2*nf2, out);
	}

	fbank_close();
	return 0;
}



