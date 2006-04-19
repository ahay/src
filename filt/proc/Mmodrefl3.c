/* Normal reflectivity modeling. 

In this version, the input contains Vp, Vs, and density into one file. 
The output contains approximate PP and PS tau-p seismograms.

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
#include <math.h>

#include <rsf.h>

#include "spline.h"

int main(int argc, char* argv[])
{
    int nt, n1, i1, i2, n2, nw, ip, np, three;
    float *a, *b, *r, *tpp, *tps, *app, *aps, *spline, **pp, **ps;
    float dt, tp,ts, a1,a2, b1,b2, r1,r2, d1, dr, da, db, ab, dp, p, as, bs, ad1, bd2; 
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&three) || three != 3) sf_error("Need n2=3 in input");
    n2 = sf_leftsize(in,2);

    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");

    if (!sf_getint("nt",&nt)) sf_error("Need nt=");
    /* time samples */
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
    /* time sampling */

    if (!sf_getint("nt",&np)) sf_error("Need np=");
    /* slope samples */
    if (!sf_getfloat("dt",&dp)) sf_error("Need dp=");
    /* slope sampling */

    if (!sf_getint("nw",&nw)) nw=4;
    /* interpolation length */

    sf_putint(out,"n1",nt);
    sf_putfloat(out,"d1",dt);
    sf_putfloat(out,"o1",0.);

    sf_putint(out,"n2",np);
    sf_putfloat(out,"d2",dp);
    sf_putfloat(out,"o2",0.);

    sf_putint(out,"n3",2);
    sf_putint(out,"n4",n2);

    a = sf_floatalloc(n1);
    b = sf_floatalloc(n1);
    r = sf_floatalloc(n1);
    tpp = sf_floatalloc(n1);
    tps = sf_floatalloc(n1);
    app = sf_floatalloc(n1);
    aps = sf_floatalloc(n1);

    spline = sf_floatalloc(nt);
    pp = sf_floatalloc2(nt,np);
    ps = sf_floatalloc2(nt,np);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(a,n1,in); /* Vp */
	sf_floatread(b,n1,in); /* Vs */
	sf_floatread(r,n1,in); /* rho */

	for (ip=0; ip < np; ip++) {
	    p = ip*dp;

	    tp = ts = 0.;
	    a2 = a[0];
	    b2 = b[0];
	    r2 = r[0];
	    for (i1=1; i1 < n1; i1++) {
		as = a2*p; if (as > 1.) sf_error("p=%g is postcritical",p);
		bs = b2*p; if (bs > 1.) sf_error("p=%g is postcritical",p);

		ad1 = d1*sqrtf(1.-as*as)/a2;
		bd1 = d1*sqrtf(1.-bs*bs)/b2;

		tp += 2.*ad1;
		ts += ad1 + bd1;
	
		a1 = a2;
		a2 = a[i1];
		b1 = b2;
		b2 = b[i1];
		r1 = r2;
		r2 = r[i1];
		
		tpp[i1] = tp;
		tps[i1] = ts;
		
		da = (a2-a1)/(a2+a1);
		db = (b2-b1)/(b2+b1);
		dr = (r2-r1)/(r2+r1);
		ab = (a2+a1)/(b2+b1);

		app[i1] = (da + dr) + (da - 4.*(2.*db+dr)/(ab*ab))*as;
		aps[i1] = 4.*db/ab + (1.+2./ab)*dr;
	    }
    
	    sf_int1_init (tpp, 0., dt, nt, sf_spline_int, nw, n1);
	    sf_int1_lop (true,false,nt,n1,spline,app);
	    spline_post(nw, 0, 1, nt, spline, trace);
	    sf_floatwrite(trace,nt,out);
	    sf_int1_lop (true,false,nt,n1,spline,bpp);
	    spline_post(nw, 0, 1, nt, spline, trace);
	    sf_floatwrite(trace,nt,out);
	    
	    sf_int1_init (tps, 0., dt, nt, sf_spline_int, nw, n1);
	    sf_int1_lop (true,false,nt,n1,spline,aps);
	    spline_post(nw, 0, 1, nt, spline, trace);
	    sf_floatwrite(trace,nt,out);
	}
    }    

    exit(0);
}

