/* Prestack shot-profile Kirchhoff migration in constant velocity. 

Requires the input to be in (time,offset,shot)
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
    bool aal, off;
    int nt,nx, nh, ix,it,ih, ns, is;
    float dt,dx, t0,x, v, dh, h0, h, t, sq, ds, s0, s, x0;
    float *time=NULL, *str=NULL, *tx=NULL, *amp=NULL, **cinp=NULL, *cout=NULL;
    sf_file inp=NULL, out=NULL;

    sf_init (argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("aal",&aal)) aal=true;
    /* if y, apply antialiasing */

    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histint(inp,"n3",&ns)) ns=1;

    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No o1= in input");
    if (!sf_histfloat(inp,"o1",&t0)) sf_error("No o1= in input");

    if (!sf_histfloat(inp,"d2",&dh)) sf_error("No d2= in input");
    if (!sf_histfloat(inp,"o2",&h0)) sf_error("No o2= in input");

    if (!sf_histfloat(inp,"d3",&ds)) sf_error("No d3= in input");
    if (!sf_histfloat(inp,"o3",&s0)) sf_error("No o3= in input");

    if (!sf_getint("nx",&nx)) nx=ns;
    if (!sf_getfloat("dx",&dx)) dx=ds;
    if (!sf_getfloat("x0",&x0)) x0=s0;

    if (!sf_getbool("offset",&off)) off=false;
    /* if y, the output is in offset */

    if (off) {
	sf_putint(out,"n3",nx);
	sf_putfloat(out,"d3",dx);
	sf_putfloat(out,"o3",x0);
	sf_putint(out,"n4",ns);
    } else {
	sf_putint(out,"n2",nx);
	sf_putfloat(out,"d2",dx);
	sf_putfloat(out,"o2",x0);
    }

    if (!sf_getfloat ("vel",&v)) sf_error("Need vel=");
    /* velocity */

    ds *= 2.0/v; s0 *= 2.0/v;
    dh *= 2.0/v; h0 *= 2.0/v;
    dx *= 2.0/v; x0 *= 2.0/v;

    time = sf_floatalloc(nt);
    str = sf_floatalloc(nt);
    tx = sf_floatalloc(nt);
    amp = sf_floatalloc(nt);

    cinp = sf_floatalloc2(nt,nh);
    cout = sf_floatalloc(nt);

    sf_aastretch_init (false, nt, t0, dt, nt);

    for (is=0; is < ns; is++) {
	sf_warning("shot %d of %d",is+1, ns);
	s = s0 + is*ds;

	sf_floatread (cinp[0],nt*nh,inp);

	for (ix = 0; ix < nx; ix++) {
	    x = x0 + ix*dx - s;

	    for (it=0; it < nt; it++) {
		t = t0 + it*dt;
		time[it] = hypotf(t,x);
		cout[it] = 0.;
	    }

	    for (ih=0; ih < nh; ih++) {
		h = x - (h0 + ih*dh);
	
		for (it=0; it < nt; it++) {
		    t = t0 + it*dt;
		    sq = hypotf(t,h);

		    str[it] = 0.5f*(time[it]+sq);
		    tx[it] = fabsf(0.5f*h*dh/sq);
		    amp[it]=1.;
		}  /* it */

		sf_aastretch_define (str, tx, amp);
		sf_aastretch_lop (true,true,nt,nt,cout,cinp[ih]);

		if (off) sf_floatwrite (cout,nt,out);
	    } /* h */

	    if (!off) sf_floatwrite (cout,nt,out);
	} /* x */
    } /* s */

    exit(0);
}
