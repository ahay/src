/* Objective function of two dips estimation with PWD filters. */
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

#include "twoplane2.h"

int main(int argc, char* argv[])
{
    bool drift;
    int i1, i2, n1, n2, np, ip, nq, iq, nw, nj, n12;
    float p0, dp, p, q0, dq, q, **pp, **qq, *xx, *yy, **obj;
    sf_file in, of;

    sf_init (argc,argv);
    in = sf_input("in");
    of = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1 in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2 in input");
    n12 = n1*n2;

    if (!sf_getint("np",&np)) np=100;
    /* number of dips */
    if (!sf_getfloat("p0",&p0)) sf_error("Need p0=");
    /* first dip origin */
    if (!sf_getfloat("dp",&dp)) dp=2*p0/(1.-np);
    /* first dip sampling */

    if (!sf_getint("nq",&nq)) nq=100;
    /* number of dips */
    if (!sf_getfloat("q0",&q0)) sf_error("Need q0=");
    /* second dip origin */
    if (!sf_getfloat("dq",&dq)) dq=2*q0/(1.-nq);
    /* second dip sampling */

    sf_putint(of,"n1",np);
    sf_putfloat(of,"d1",dp);
    sf_putfloat(of,"o1",p0);

    sf_putint(of,"n2",nq);
    sf_putfloat(of,"d2",dq);
    sf_putfloat(of,"o2",q0);

    if (!sf_getint("order",&nw)) nw=1;
    /* [1,2,3] accuracy order */
    if (nw < 1 || nw > 3) 
	sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);
    if (!sf_getint("nj",&nj)) nj=1;
    /* antialiasing */

    if (!sf_getbool("drift",&drift)) drift=false;
    /* if shift filter */

    xx = sf_floatalloc(n12);
    yy = sf_floatalloc(n12);
    pp = sf_floatalloc2(n1,n2);
    qq = sf_floatalloc2(n1,n2);
    obj = sf_floatalloc2(np,nq);

    sf_floatread (xx,n12,in);
    twoplane2_init (nw,nj,nj,n1,n2,drift,pp,qq);

    for (iq=0; iq < nq; iq++) {
        q = q0 + iq*dq;
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		qq[i2][i1] = q;
	    }
	}

	for (ip=0; ip < np; ip++) {
	    p = p0 + ip*dp;
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    pp[i2][i1] = p;
		}
	    }

	    twoplane2_lop(false,false,n12,n12,xx,yy);

	    obj[iq][ip] = 0.;
	    for (i1=0; i1 < n12; i1++) {
		obj[iq][ip] += yy[i1]*yy[i1];
	    }
	}
    }

    sf_floatwrite(obj[0],np*nq,of);

    exit(0);
}
