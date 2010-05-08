/* Hyperbolic Radon transform */
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
#include "veltran.h"

int main(int argc, char* argv[])
{
    int nt, n3, nx, nv, i3, ntx, ntv;
    bool adj, pull;
    float o1,d1, x0,dx, v0,dv, anti, s1,s0,ds;
    float *cmp=NULL, *vscan=NULL;
    sf_file in=NULL, out=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");

    n3 = sf_leftsize(in,2);

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */
    if (!sf_getfloat("anti",&anti)) anti=1.;
    /* antialiasing */

    if (adj) {
	if (!sf_histint(in,"n2",&nx))   sf_error("No n2= in input");
	if (!sf_histfloat(in,"o2",&x0)) sf_error("No o2= in input");
	if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");

	if (!sf_getint("nv",&nv) && !sf_histint(in,"nv",&nv)) nv=nx;
	if (!sf_getfloat("v0",&v0) && !sf_histfloat(in,"v0",&v0))
	    sf_error("need v0=");
	if (!sf_getfloat("dv",&dv) && !sf_histfloat(in,"dv",&dv))
	    sf_error("need dv=");

	sf_putfloat(out,"x0",x0);
	sf_putfloat(out,"dx",dx);
	sf_putint(out,"nx",nx);

	sf_putfloat(out,"v0",v0);
	sf_putfloat(out,"dv",dv);
    } else {
	if (!sf_histint(in,"n2",&nv))   sf_error("No n2= in input");
	if (!sf_histfloat(in,"v0",&v0)) sf_error("No v0= in input");
	if (!sf_histfloat(in,"dv",&dv)) sf_error("No dv= in input");

	if (!sf_getint("nx",&nx) && !sf_histint(in,"nx",&nx)) nx=nv;
	if (!sf_getfloat("x0",&x0) && !sf_histfloat(in,"x0",&x0))
	    sf_error("need x0=");
	if (!sf_getfloat("dx",&dx) && !sf_histfloat(in,"dx",&dx))
	    sf_error("need dx=");
    }

    if (!sf_getfloat("s02",&s1)) s1=0.;
    /* reference slowness squared (for antialiasing) */

    s0 = 1./(v0 +(nv-1)*dv); s0 *= s0; ds = (1./(v0*v0) - s0)/(nv-1);
/*    s0 = ds;		               ds = (1./(v0*v0) - s0)/(nv-1); */

    if (adj) {
	sf_putstring(out,"label2","Slowness Squared");
	sf_putstring(out,"unit2","s\\^2\\_/km\\^2");
	sf_putfloat(out,"o2",s0);
	sf_putfloat(out,"d2",ds);
	sf_putint(out,"n2",nv);
    } else {
	sf_putstring(out,"label2","Offset");
	sf_putstring(out,"unit2","km");
	sf_putfloat(out,"o2",x0);
	sf_putfloat(out,"d2",dx);
	sf_putint(out,"n2",nx);
    }

    ntx = nt*nx;
    ntv = nt*nv;

    cmp = sf_floatalloc(ntx);
    vscan = sf_floatalloc(ntv);

    if (!sf_getbool("pull",&pull)) pull=true;
    /* pull or push operator */

    veltran_init (pull, x0, dx, nx, s0, ds, nv, o1, d1, nt, 
		  s1, anti, 0, 0);

    for (i3=0; i3 < n3; i3++) { 
	if( adj) {
	    sf_floatread(cmp,ntx,in);
	} else {
	    sf_floatread(vscan,ntv,in);
	} 

	veltran_lop(adj,false,ntv,ntx,vscan,cmp);

	if( adj) {
	    sf_floatwrite(vscan,ntv,out);
	} else {
	    sf_floatwrite(cmp,ntx,out);
	}
    }

    exit(0);
}
