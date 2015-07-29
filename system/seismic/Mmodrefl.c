/* Normal reflectivity modeling. */
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
    int nt, n1, i1, nw;
    float *z=NULL, *a=NULL, *b=NULL, *r=NULL, *tpp=NULL, *tps=NULL, *app=NULL, *aps=NULL, *spline=NULL, *trace=NULL;
    float dt, tp,ts, a1,a2, b1,b2, r1,r2;
    sf_file depth=NULL, vp=NULL, vs=NULL, rho=NULL, dat=NULL;

    sf_init(argc,argv);

    depth = sf_input("in");
    vp = sf_input("vp");
    vs = sf_input("vs");
    rho = sf_input("rho");
    dat = sf_output("out");

    if (!sf_histint(depth,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt=");
    /* time samples */
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
    /* time sampling */
    if (!sf_getint("nw",&nw)) nw=4;
    /* interpolation length */


    sf_putint(dat,"n1",nt);
    sf_putfloat(dat,"d1",dt);
    sf_putfloat(dat,"o1",0.);

    sf_putint(dat,"n2",2);

    z = sf_floatalloc(n1);
    a = sf_floatalloc(n1+1);
    b = sf_floatalloc(n1+1);
    r = sf_floatalloc(n1+1);
    tpp = sf_floatalloc(n1);
    tps = sf_floatalloc(n1);
    app = sf_floatalloc(n1);
    aps = sf_floatalloc(n1);

    sf_floatread(z,n1,depth);
    sf_floatread(a,n1+1,vp);
    sf_floatread(b,n1+1,vs);
    sf_floatread(r,n1+1,rho);

    tp = ts = 0.;
    a2 = a[0];
    b2 = b[0];
    r2 = r[0];
    for (i1=0; i1 < n1; i1++) {
	tp += 2.*z[i1]/a2;
	ts += z[i1]/a2 + z[i1]/b2;
	
	a1 = a2;
	a2 = a[i1+1];
	b1 = b2;
	b2 = b[i1+1];
	r1 = r2;
	r2 = r[i1+1];
	
	tpp[i1] = tp;
	tps[i1] = ts;
	app[i1] = (a2-a1)/(a2+a1) + (r2-r1)/(r2+r1);
	aps[i1] = 4.*(b2-b1)/(a2+a1) + 
	    (1.+2.*(b2+b1)/(a2+a1))*(r2-r1)/(r2+r1);
    }

    spline = sf_floatalloc(nt);
    trace = sf_floatalloc(nt);

    sf_int1_init (tpp, 0., dt, nt, sf_spline_int, nw, n1);
    sf_int1_lop (true,false,nt,n1,spline,app);
    sf_spline_post(nw, 0, 1, nt, spline, trace);
    sf_floatwrite(trace,nt,dat);

    sf_int1_init (tps, 0., dt, nt, sf_spline_int, nw, n1);
    sf_int1_lop (true,false,nt,n1,spline,aps);
    sf_spline_post(nw, 0, 1, nt, spline, trace);
    sf_floatwrite(trace,nt,dat);

    exit(0);
}
