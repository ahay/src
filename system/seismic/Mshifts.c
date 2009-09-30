/* Multiple shifts. */
/*
  Copyright (C) 2009 University of Texas at Austin

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
#include "fint1.h"

int main(int argc, char* argv[])
{
    int i1, n1, i2, n2, i3, n3, ip, np, it;
    float **dat=NULL, *trace=NULL, dp, p0, p, t;
    vint1 str;
    sf_file inp=NULL, out=NULL;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&n2) || n2 < 2) sf_error("No n2= in input");
    n3 = sf_leftsize(inp,2);

    if (!sf_getint("np",&np)) sf_error("Need np=");
    /* number of slopes */
    if (!sf_getfloat("dp",&dp)) sf_error("Need dp=");
    /* slope sampling */
    if (!sf_getfloat("p0",&p0)) sf_error("Need p0=");
    /* first slope */

    sf_putint(out,"n3",np);
    sf_putfloat(out,"d3",dp);
    sf_putfloat(out,"o3",p0);

    sf_shiftdim(inp,out,3);

    dat = sf_floatalloc2(n1,n2);
    trace = sf_floatalloc(n2);
    trace[n2-1] = 0.;

    str = vint1_init (4,n1,n2-1);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(dat[0],n1*n2,inp);
	vint1_set(str,dat+1);

	for (ip=0; ip < np; ip++) {
	    p = p0 + ip*dp;

	    for (i1=0; i1 < n1; i1++) {
		t = i1+p;
		it = floorf(t);

		if (it < -1 || it > n1) {
		    for (i2=0; i2 < n2-1; i2++) {
			trace[i2]=0.;
		    }
		} else {
		    vint1_apply(str,it,t-it,false,trace);
		}

		for (i2=0; i2 < n2; i2++) {
		    dat[i2][i1] = trace[i2];
		}
	    }
	    sf_floatwrite(dat[0],n1*n2,out);
	}
    }

    exit(0);
}
