/* Normal reflectivity modeling. 

In this version, the input contains Vp, Vs, and density into one file. 
The output contains PP and PS tau-p seismograms.

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

#include "zoeppritz.h"

int main(int argc, char* argv[])
{
    bool moveout;
    int nt, n1, i1, i2, n2, ns, n1s, ip, np, three;
    float *a=NULL, *b=NULL, *r=NULL, *tpp=NULL, *tps=NULL, *app=NULL, *aps=NULL, **pp=NULL, **ps=NULL;
    float dt, tp,ts, a1,a2, b1,b2, r1,r2, eps, rc[4], ang[4];
    float d1, p0, dp, p, as, bs, ad1, bd1; 
    sf_map4 map;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&three) || three != 3)
	sf_error("Need n2=3 in input");
    n2 = sf_leftsize(in,2);

    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");

    if (!sf_getint("sparse",&ns)) ns=10;
    /* sparseness of reflectivity */

    if (!sf_getbool("moveout",&moveout)) moveout=true;
    /* if apply moveout */

    if (!sf_getint("nt",&nt)) sf_error("Need nt=");
    /* time samples */
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
    /* time sampling */

    if (!sf_getint("np",&np)) sf_error("Need np=");
    /* slope samples */
    if (!sf_getfloat("dp",&dp)) sf_error("Need dp=");
    /* slope sampling */
    if (!sf_getfloat("p0",&p0)) sf_error("Need p0=");
    /* slope origin */

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    sf_putint(out,"n1",nt);
    sf_putfloat(out,"d1",dt);
    sf_putfloat(out,"o1",0.);
    sf_putstring(out,"label1","Time");
    sf_putstring(out,"unit1","s");

    sf_putint(out,"n2",np);
    sf_putfloat(out,"d2",dp);
    sf_putfloat(out,"o2",p0);

    sf_putint(out,"n3",2);
    sf_putint(out,"n4",n2);    

    a = sf_floatalloc(n1);
    b = sf_floatalloc(n1);
    r = sf_floatalloc(n1);

    n1s = (n1-1)*ns;
    d1 /= ns;

    tpp = sf_floatalloc(n1s);
    tps = sf_floatalloc(n1s);
    app = sf_floatalloc(n1s);
    aps = sf_floatalloc(n1s);

    pp = sf_floatalloc2(nt,np);
    ps = sf_floatalloc2(nt,np);

    map = sf_stretch4_init (nt, 0., dt, n1s, eps);

    for (i2=0; i2 < n2; i2++) {
	sf_warning("CMP %d of %d",i2+1,n2);

	sf_floatread(a,n1,in); /* Vp */
	sf_floatread(b,n1,in); /* Vs */
	sf_floatread(r,n1,in); /* rho */

	for (ip=0; ip < np; ip++) {
	    p = p0+ip*dp;

	    tp = ts = 0.;
	    a2 = a[0];
	    b2 = b[0];
	    r2 = r[0];
	    for (i1=1; i1 <= n1s; i1++) {
		as = a2*p;
		if (fabsf(as) > 1.) 
		    sf_error("p=%g is postcritical (vp=%g)",p,a2);

		bs = b2*p; 
		if (fabsf(bs) > 1.) 
		    sf_error("p=%g is postcritical (vs=%g)",p,b2);
		
		ad1 = d1*sqrtf(1.-as*as)/a2;
		bd1 = d1*sqrtf(1.-bs*bs)/b2;

		if (moveout) {
		    tp += 2.*ad1;
		    ts += ad1 + bd1;
		} else {
		    tp += 2.*d1/a2;
		    ts += d1/a2+d1/b2;
		}

		tpp[i1-1] = tp;
		tps[i1-1] = ts;
	
		a1 = a2;
		b1 = b2;
		r1 = r2;

		if (0==i1%ns) {
		    a2 = a[i1/ns];
		    b2 = b[i1/ns];
		    r2 = r[i1/ns];
		    zoeppritz (4,a1,a2,b1,b2,r1,r2,true,p,rc,ang);
		} else {
		    rc[0] = rc[1] = 0.;
		}

		app[i1-1] = rc[0];
		aps[i1-1] = rc[1];
	    }

	    sf_stretch4_define (map,tpp);
	    sf_stretch4_apply (false,map,app,pp[ip]);

	    sf_stretch4_define (map,tps);
	    sf_stretch4_apply (false,map,aps,ps[ip]);
	}

	sf_floatwrite(pp[0],nt*np,out);
	sf_floatwrite(ps[0],nt*np,out);
    }

    exit(0);
}
