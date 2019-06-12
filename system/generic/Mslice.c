/* Extract a slice using picked surface (usually from a stack or a semblance).

See also: sfpick.

June 2019 program of the month:
http://ahay.org/blog/2019/06/12/program-of-the-month-sfslice/
*/
/*
  Copyright (C) 2004 University of Texas at Austin

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.h>

int main(int argc, char* argv[])
{
    int nt, ns, nx, ix, is, it;
    float *trace=NULL, *picks=NULL, **semb=NULL, s0, ds, s;
    sf_file in=NULL, out=NULL, pick=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    pick = sf_input("pick");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");

    if (!sf_histint(in,"n2",&ns)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"o2",&s0)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"d2",&ds)) sf_error("No d2= in input");

    nx = sf_unshiftdim(in,out,2);

    semb = sf_floatalloc2(nt,ns);
    picks = sf_floatalloc (nt);
    trace = sf_floatalloc (nt);

    for (ix=0; ix < nx; ix++) {
	sf_floatread (semb[0],nt*ns,in);
	sf_floatread (picks,nt,pick);

	for (it=0; it < nt; it++) {
	    s = (picks[it] - s0)/ds;
	    is = floorf(s); s -= is;
	    if (is >= 0 && is < ns-1) {
		trace[it] = s * semb[is+1][it] + (1.-s) * semb[is][it];
	    } else {
		trace[it] = 0.;
	    }
	}

	sf_floatwrite (trace,nt,out);
    }

    exit (0);
}
