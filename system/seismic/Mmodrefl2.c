/* Normal reflectivity modeling. 

In this version, the input contains Vp, Vs, and density into one file. 
The output contains PP intercept, PP gradient, and PS gradient.

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

int main(int argc, char* argv[])
{
    int nt, n1, i1, i2, n2, nw, three;
    float *a=NULL, *b=NULL, *r=NULL, *tpp=NULL, *tps=NULL, *app=NULL, *aps=NULL, *bpp=NULL, *spline=NULL, *trace=NULL;
    float dt, tp,ts, a1,a2, b1,b2, r1,r2, d1, dr, da, db, ab; 
    sf_file in=NULL, out=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&three) || three != 3)
	sf_error("Need n2=3 in input");
    n2 = sf_leftsize(in,2);

    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");

    if (!sf_getint("nt",&nt)) sf_error("Need nt=");
    /* time samples */
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
    /* time sampling */
    if (!sf_getint("nw",&nw)) nw=4;
    /* interpolation length */

    sf_putint(out,"n1",nt);
    sf_putfloat(out,"d1",dt);
    sf_putfloat(out,"o1",0.);

    a = sf_floatalloc(n1);
    b = sf_floatalloc(n1);
    r = sf_floatalloc(n1);
    tpp = sf_floatalloc(n1);
    tps = sf_floatalloc(n1);
    app = sf_floatalloc(n1);
    bpp = sf_floatalloc(n1);
    aps = sf_floatalloc(n1);

    spline = sf_floatalloc(nt);
    trace = sf_floatalloc(nt);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(a,n1,in);
	sf_floatread(b,n1,in);
	sf_floatread(r,n1,in);

	tp = ts = 0.;
	a2 = a[0];
	b2 = b[0];
	r2 = r[0];
	for (i1=1; i1 < n1; i1++) {
	    tp += 2.*d1/a2;
	    ts += d1/a2 + d1/b2;
	
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

	    app[i1] = da + dr;
	    bpp[i1] = da - 4.*(2.*db+dr)/(ab*ab);
	    aps[i1] = 4.*db/ab + (1.+2./ab)*dr;
	}

	sf_int1_init (tpp, 0., dt, nt, sf_spline_int, nw, n1);
	sf_int1_lop (true,false,nt,n1,spline,app);
	sf_spline_post(nw, 0, 1, nt, spline, trace);
	sf_floatwrite(trace,nt,out);
	sf_int1_lop (true,false,nt,n1,spline,bpp);
	sf_spline_post(nw, 0, 1, nt, spline, trace);
	sf_floatwrite(trace,nt,out);

	sf_int1_init (tps, 0., dt, nt, sf_spline_int, nw, n1);
	sf_int1_lop (true,false,nt,n1,spline,aps);
	sf_spline_post(nw, 0, 1, nt, spline, trace);
	sf_floatwrite(trace,nt,out);
    }

    exit(0);
}
