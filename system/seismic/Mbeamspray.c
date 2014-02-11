/* 2-D beam spraying. */
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

#include "aastretch2.h"

int main(int argc, char* argv[])
{
    int n1, nc, nd, n3, i3, nb, id, ic, i1, ib;
    float **dense, *a, *p, *c, *time, *delt1, *delt2, *ampl, d2, f, eps;
    sf_file in, out, dip, cur;
    float f1=6.77917496069906;
    float f2=1.75843756426359;

    sf_init(argc,argv);
    in = sf_input("in");
    dip = sf_input("dip");
    cur = sf_input("cur");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nc)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_histfloat(in,"d2",&d2)) d2=1.;

    if (!sf_getint("rect",&nb)) nb=3;
    /* smoothing radius */

    if (!sf_getfloat("eps",&eps)) eps=1.0;
    /* experimental */

    nd = (nc-1)*nb+1;
    sf_putint(out,"n2",nd);
    sf_putfloat(out,"d2",d2/nb);

    dense = sf_floatalloc2(n1,nd);
    a = sf_floatalloc(n1);
    p = sf_floatalloc(n1);
    c = sf_floatalloc(n1);

    time = sf_floatalloc(n1);
    ampl = sf_floatalloc(n1);
    delt1 = sf_floatalloc(n1);
    delt2 = sf_floatalloc(n1);

    /*** Initialize stretch ***/
    aastretch2_init (n1, 0., 1.0, n1);

    f = 0.25*(1+0.25*SF_PI*SF_PI);
    f *= f/(nb*nb);
    f = 9/(8*SF_PI*eps)/(nb*nb);
    f1 *= f;
    f2 *= f;

    for (i3=0; i3 < n3; i3++) {
	for (id=0; id < nd; id++) {
	    for (i1=0; i1 < n1; i1++) {
		dense[id][i1] = 0.;
	    }
	}
	for (ic=0; ic < nc; ic++) {
	    sf_floatread(a,n1,in);
	    sf_floatread(p,n1,dip);
	    sf_floatread(c,n1,cur);

	    for (ib=-2*nb; ib <= 2*nb; ib++) {
		id=ic*nb+ib;
		
		if (id >= nd) break;
		if (id < 0) continue;
		
		for (i1=0; i1 < n1; i1++) {
		    time[i1] = i1+(p[i1]+c[i1]*ib/2)*ib;
		    delt1[i1] = f1*ib*ib+1.0;
		    delt2[i1] = f2*ib*ib+1.0;
		    ampl[i1] = 0.4;
		}
		
		aastretch2_define (time,delt1,delt2,ampl);
		aastretch2_lop (false,true,n1,n1,a,dense[id]);
	    }
	}
	sf_floatwrite(dense[0],n1*nd,out);
    }

    exit(0);
}

