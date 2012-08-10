/* linear phase filter */

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

#include "lphase.h"

int main(int argc, char* argv[])
{
	sf_file out;
	int nf, interp;
	float **c;

	sf_init(argc, argv);
	out = sf_output ("out");

	if(!sf_getint("order", &nf)) nf=1;
	/* order of linear phase filter */
	if (!sf_getint("interp",&interp)) interp=0;
	/* interpolation method: 
	0: maxflat
	1: Lagrange 
	2: B-Spline */

	sf_putint(out, "n1", 2*nf+1);
	sf_putfloat(out, "o1", -nf);
	sf_putfloat(out, "d1", 1);
	sf_putint(out, "n2", 2*nf+1);
	sf_putfloat(out, "o2", 0);
	sf_putfloat(out, "d2", 1);

	c = lphase(nf, interp);
	sf_floatwrite(c[0], (2*nf+1)*(2*nf+1), out);

	free(c[0]);
	free(c);
	return 0;
}

