/* Objective function of dip estimation with PWD filters. */
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

#include "allp2.h"

int main(int argc, char* argv[])
{
    bool drift;
    int i1, i2, n1, n2, np, ip, nw, nj;
    float p0, dp, p, **pp, **xx, **yy, *obj;
    allpas2 ap;
    sf_file in, of;

    sf_init (argc,argv);
    in = sf_input("in");
    of = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1 in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2 in input");

    if (!sf_getint("np",&np)) np=100;
    /* number of dips */
    if (!sf_getfloat("p0",&p0)) sf_error("Need p0=");
    /* dip origin */
    if (!sf_getfloat("dp",&dp)) dp=2*p0/(1.-np);
    /* dip sampling */

    if (!sf_getbool("drift",&drift)) drift=false;
    /* if shift filter */

    sf_putint(of,"n1",np);
    sf_putfloat(of,"d1",dp);
    sf_putfloat(of,"o1",p0);
    sf_putint(of,"n2",1);

    if (!sf_getint("order",&nw)) nw=1;
    /* [1,2,3] accuracy order */
    if (nw < 1 || nw > 3) 
	sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);
    if (!sf_getint("nj",&nj)) nj=1;
    /* antialiasing */

    xx = sf_floatalloc2(n1,n2);
    yy = sf_floatalloc2(n1,n2);
    pp = sf_floatalloc2(n1,n2);
    obj = sf_floatalloc(np);

    sf_floatread (xx[0],n1*n2,in);
    ap = allpass2_init (nw,nj,n1,n2,drift,pp);

    for (ip=0; ip < np; ip++) {
        p = p0 + ip*dp;
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		pp[i2][i1] = p;
	    }
	}

	allpass21 (false, ap,xx,yy);

        obj[ip] = 0.;
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		obj[ip] += yy[i2][i1]*yy[i2][i1];
	    }
	}
    }

    sf_floatwrite(obj,np,of);

    exit(0);
}
