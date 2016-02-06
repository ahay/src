/* Velocity-independent phase-space zero-offset migration. */
/*
  Copyright (C) 2015 University of Texas at Austin

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

int main(int argc, char* argv[])
{
    int nt, nx, np, it, ix, ip, ixin, ix1, ix2, n12;
    float t0, dt, dx, p0, dp, t, x, p, px, vmin, vmax;
    float **img, *amp, *tx, *str, *add, **cinp, **cout;
    sf_file inp, mig;

    sf_init (argc,argv);
    inp = sf_input("in");
    mig = sf_output("out");

    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(inp,"n3",&np)) sf_error("No n3= in input");

    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(inp,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histfloat(inp,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(inp,"o3",&p0)) sf_error("No o3= in input");
    if (!sf_histfloat(inp,"d3",&dp)) sf_error("No d3= in input");

    if (!sf_getfloat("vmin",&vmin)) sf_error("Need vmin=");
    if (!sf_getfloat("vmax",&vmax)) sf_error("Need vmax=");

    /* convert to slowness */
    vmin = 1.0f/vmin;
    vmax = 1.0f/vmax;

    img = sf_floatalloc2(nt,nx);
    amp = sf_floatalloc(nt);
    tx = sf_floatalloc(nt);
    str = sf_floatalloc(nt);
    add = sf_floatalloc(nt);
    
    cinp = sf_floatalloc2(nt,nx);
    cout = sf_floatalloc2(nt,nx);

    n12 = nt*nx;

    sf_aastretch_init (false, nt, t0, dt, nt);

    for (ip=0; ip < np; ip++) {
	sf_warning("slope %d of %d;",ip+1,np);

	p = p0+ip*dp;
	sf_floatread (cinp[0],n12,inp);
	
	for (it=0; it<n12; it++) {
	    cout[0][it] = 0.0f;
	}

	if (fabsf(p) < vmax) {
	
	    /* loop over shifts */
	    for (ix = 1-nx; ix <= nx-1; ix++) {
		x = ix*dx;

		px = p*x;

		if (ix > 0) {
		    ix1 = 0;
		    ix2 = nx-ix;
		} else {
		    ix1 = -ix;
		    ix2 = nx;
		}
		
		for (it=0; it < nt; it++) {
		    t = t0+it*dt;

		    str[it] = (px+hypotf(2*t,px))/2;

		    /* antialiasing */
		    tx[it]=fabsf(p*dx*(1+px/hypotf(2*t,px))/2);
		    
		    if ((p > vmin*vmin*x/hypotf(t,vmin*x) && p < vmax*vmax*x/hypotf(t,vmax*x)) ||
			(p < vmin*vmin*x/hypotf(t,vmin*x) && p > vmax*vmax*x/hypotf(t,vmax*x))) {
			amp[it] = 1.0f;
		    } else {
			amp[it] = 0.0f;
		    }
		}
		
		sf_aastretch_define (str, tx, amp);
		
		for (ixin=ix1; ixin < ix2; ixin++) {
		    sf_aastretch_lop (false,false,nt,nt, cinp[ixin], add);
		    for (it=0; it < nt; it++) {
			cout[ixin+ix][it] += add[it];
		    }
		}

	    }
	}

	sf_floatwrite (cout[0],n12,mig);
    }
    sf_warning(".");

    exit(0);
}
	
