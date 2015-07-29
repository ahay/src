/* Angle-gather constant-velocity time migration. */
/*
  Copyright (C) 2006 University of Texas at Austin

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

int main(int argc, char *argv[])
{
    int nt,nx,nh,na,ng, ix,iy,ik,iz,it,ih,ia,ig;
    float dt,dx,t0,x, vel, t,z, tf,hf,xf,den;
    float dh,h0,h, g0,dg,g, a,da, ca,sa,cg,sg, wt, amax;
    float ***dat=NULL, **img=NULL;
    sf_file in=NULL, out=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&nh)) sf_error("No n3= in input");

    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"d3",&dh)) sf_error("No d3= in input");

    if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"o3",&h0)) sf_error("No o3= in input");

    if (!sf_getfloat("vel",&vel)) sf_error("Need vel=");
    /* velocity */

    if (!sf_getint("ng",&ng)) sf_error("Need ng=");
    /* number of reflection angles */

    if (!sf_getfloat("dg",&dg)) sf_error("Need dg=");
    /* reflection angle sampling */

    if (!sf_getfloat("g0",&g0)) sf_error("Need g0=");
    /* reflection angle origin */

    sf_putint(out,"n3",ng);
    sf_putfloat(out,"d3",dg);
    sf_putfloat(out,"o3",g0);

    if (!sf_getint("na",&na)) na=nx;
    /* number of dip angles */
    if (!sf_getfloat("a",&amax)) amax=80.;
    /* maximum dip angle */

    /* angles to radians */
    dg   *= SF_PI / 180.;
    g0   *= SF_PI / 180.;
    amax *= SF_PI / 180.;

    dat = sf_floatalloc3(nt,nx,nh);
    img = sf_floatalloc2(nt,nx);

    sf_floatread (dat[0][0],nt*nx*nh,in);

    for (ig=0; ig < ng; ig++) {
	sf_warning("angle %d of %d;",ig+1,ng);
	g = g0 + ig * dg;
	for (ix=0; ix < nx; ix++) {
	    for (it=0; it < nt; it++) {
		img[ix][it] = 0.;
	    }
	}
	if (fabsf(g) >= amax) {
	    sf_floatwrite(img[0],nt*nx,out);
	    continue;
	}
	da = (amax - fabsf(g))/na;
	cg = cosf(g);
	sg = sinf(g);

	for (ia=0; ia < na; ia++) {
	    a = ia*da;
	    ca = cosf(a);
	    sa = sinf(a);

	    den = ca*ca - sg*sg;
	    tf = ca*cg/den;
	    hf = (vel*0.5) * sg*cg/den;
	    xf = (vel*0.5) * sa*ca/den;

	    for (iz=0; iz < nt; iz++) {
		z = t0 + iz*dt;
		t = (z*tf-t0)/dt; it = floorf(t); t -= it;
		if (it >= nt-1) break;
		if (it < 0) continue;

		h = (z*hf-h0)/dh; ih = floorf(h); h -= ih; 
		if (ih >= nh-1) break;
		if (ih < 0) continue;

		x = z*xf/dx; ix = floorf(x); x -= ix;

		wt = (da * dg) * (z / sqrtf(dx * dh)) * sa *
		    (ca * ca + sg *sg)/(den*den*sqrtf(den));
  
		for(iy=0; iy < nx-ix-1; iy++) {
		    ik = iy + ix;
		    img[iy][iz]
			+= (1.-t)*(1.-x)*(1.-h)*wt * dat[ih  ][ik  ][it  ] 
			+      t *(1.-x)*(1.-h)*wt * dat[ih  ][ik  ][it+1]
			+  (1.-t)*    x *(1.-h)*wt * dat[ih  ][ik+1][it  ]
			+      t *    x *(1.-h)*wt * dat[ih  ][ik+1][it+1]
			+  (1.-t)*(1.-x)*    h *wt * dat[ih+1][ik  ][it  ]
			+      t *(1.-x)*    h *wt * dat[ih+1][ik  ][it+1]
			+  (1.-t)*    x *    h *wt * dat[ih+1][ik+1][it  ]
			+      t *    x *    h *wt * dat[ih+1][ik+1][it+1];
		}
		for (iy=ix+1; iy < nx; iy++) {
		    ik = iy - ix - 1;
		    img[iy][iz]
		    += (1.-t)*(1.-x)*(1.-h)*wt * dat[ih  ][ik  ][it  ] 
			+      t *(1.-x)*(1.-h)*wt * dat[ih  ][ik  ][it+1]
			+  (1.-t)*    x *(1.-h)*wt * dat[ih  ][ik+1][it  ]
			+      t *    x *(1.-h)*wt * dat[ih  ][ik+1][it+1]
			+  (1.-t)*(1.-x)*    h *wt * dat[ih+1][ik  ][it  ]
			+      t *(1.-x)*    h *wt * dat[ih+1][ik  ][it+1]
			+  (1.-t)*    x *    h *wt * dat[ih+1][ik+1][it  ]
			+      t *    x *    h *wt * dat[ih+1][ik+1][it+1];
		}
	    } /* iz */
	} /* ia */

	sf_floatwrite(img[0],nt*nx,out);
    } /* ig */
    sf_warning(".");

    exit(0);
}
