/* Linear PHase filter COEFficients */

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
#include "lphpoly.h"

int main(int argc, char* argv[])
{
	sf_file out;
	int nf, m, n;
	float **c;
	char *interp;

	sf_init(argc, argv);
	out = sf_output ("out");

	if(!sf_getint("m", &m)) m=1;
	/* b[-m, ... ,n] */
	if(!sf_getint("n", &n)) n=1;
	/* b[-m, ... ,n] */
	if ((interp=sf_getstring("interp"))==NULL) interp="maxflat";
	/* interpolation method: maxflat lagrange bspline */

	nf = m+n+1;
	sf_putint(out, "n1", nf);
	sf_putint(out, "o1", -m);
	sf_putint(out, "d1", 1);
	sf_putint(out, "n2", nf);
	sf_putint(out, "o2", 0);
	sf_putint(out, "d2", 1);

	c = lphpoly(m, n, interp);
	if(c==NULL) sf_error("interp=%s incorrect", interp);
	sf_floatwrite(c[0], nf*nf, out);

	free(c[0]);
	free(c);
	return 0;
}

