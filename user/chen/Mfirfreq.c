/* frequency response of fir filters */

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

#include "dsp.h"

int main(int argc, char* argv[])
{
	sf_file in, out;
	int n1, n2, mw;
	int iw, i2;
	float **c, ow, dw, o1, d1;
	sf_complex **d;

	sf_init(argc, argv);
	in  = sf_input ("in");
	out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
	if (!sf_histfloat(in,"o1", &o1)) sf_error("o1");
	if (!sf_histfloat(in,"d1", &d1)) sf_error("d1");

	if(!sf_getint("nw", &mw)) mw=101;
	/* samples in frequency domain */
	if(!sf_getfloat("ow", &ow)) ow=-0.5;
	/* first frequency */
	if(!sf_getfloat("dw", &dw)) dw=0.01;
	/* frequency increment*/

	sf_putint(out, "n1", mw);
	sf_putfloat(out, "o1", ow);
	sf_putfloat(out, "d1", dw);
	sf_settype(out,SF_COMPLEX);

	c = sf_floatalloc2(n1, n2);
	d = sf_complexalloc2(mw, n2);

	sf_floatread(c[0], n1*n2, in);

	for(i2=0; i2<n2; i2++)
	for(iw=0; iw<mw; iw++)
	d[i2][iw] = fir_freq((int)(o1-0.5), (int)(o1+(n1-1)+0.5), 
		c[i2]-(int)(o1-0.5), (ow+dw*iw));

	sf_complexwrite(d[0], mw*n2, out);

	free(c[0]);
	free(c);
	free(d[0]);
	free(d);
	return 0;
}

