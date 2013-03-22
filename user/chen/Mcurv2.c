/* Joint estimation of curvature and slope */

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

#include "curv2.h"

int main(int argc,char*argv[])
{
	sf_file in, out, slope;
	float **u, **zx;
	int n1, n2, n3, i3;
	int m, n, niter;
	char *interp;
	sf_init(argc,argv);

	in  = sf_input("in");
	out = sf_output("out");

	slope = sf_output("slope");

	if (!sf_histint(in, "n1", &n1)) sf_error("n1 needed");
	if (!sf_histint(in, "n2", &n2)) sf_error("n2 needed");
	if (!sf_histint(in, "n3", &n3)) n3=1;
	if (!sf_getint("niter", &niter)) niter=5;
	/* iterations */

    if(!sf_getint("m", &m)) m=1;
    /* b[-m, ... ,n] */
    if(!sf_getint("n", &n)) n=1;
    /* b[-m, ... ,n] */
    if ((interp=sf_getstring("interp"))==NULL) interp="maxflat";
    /* interpolation method: maxflat lagrange bspline */

	u  = sf_floatalloc2(n1, n2);
	zx  = sf_floatalloc2(n1, n2);

	curv2_init(m, n, interp, n1, n2, niter);

	for(i3=0; i3<n3; i3++)
	{
		sf_floatread(u[0], n1*n2, in);
		curv2(u, zx);
		sf_floatwrite(u[0], n1*n2, out);
		sf_floatwrite(zx[0], n1*n2, slope);
	}
	curv2_close();
	free(*u);
	free(u);
	free(*zx);
	free(zx);
	exit(0);
}

