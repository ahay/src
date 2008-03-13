/* Plane-wave destruction with factored helix filters. */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

#include "nhelix.h"
#include "nhelicon.h"
#include "wilson.h"
#include "regrid.h"
#include "compress.h"

int main(int argc, char* argv[]) 
{
    int ntxy, nt, nx, ny, pt, px, i, niter, npx, npy, np, ix, iy, ptx, ia, m[3], n[3];
    float a0, s0, p, q, f, dt, dx, dy, pmin, pmax, qmin, qmax, dp, dq, eps, *pp, *qq; 
    sf_filter ss, aa;
    nfilter pfilt;
    sf_file in, out, dip;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    dip = sf_input("dip");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&ny)) sf_error("No n3= in input");
    ntxy = nt*nx*ny;
    n[0] = nt; n[1] = nx; n[2] = ny;

    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"d3",&dy)) sf_error("No d3= in input");

    if (!sf_getfloat("eps",&eps)) eps=0.001;
    if (!sf_getint("nt",&pt));
    if (!sf_getint("nx",&px));
    ptx = pt*px;
    m[0] = pt; m[1] = px; m[2] = px;

    if (!sf_getint("npx",&npx)) npx=100;
    if (!sf_getint("npy",&npy)) npy=100;
    np = npx *npy;

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of iterations */

    ss = sf_allocatehelix(8);
    ss->lag[0] = 1;
    ss->lag[1] = 2;
    ss->lag[2] = pt - 1;
    ss->lag[3] = pt;
    ss->lag[4] = pt + 1;
    ss->lag[5] = pt * px - 1;
    ss->lag[6] = pt * px;
    ss->lag[7] = pt * px + 1;

    pp = sf_floatalloc(ntxy);
    qq = sf_floatalloc(ntxy);

    sf_floatread(pp,ntxy,dip);
    sf_floatread(qq,ntxy,dip);
    sf_fileclose(dip);

    pfilt = (nfilter) sf_alloc(1,sizeof(*pfilt));
    pfilt->hlx = (sf_filter*) sf_alloc(np,sizeof(sf_filter));
    pfilt->pch = sf_intalloc(ntxy);
    pfilt->mis = NULL;

    pmin = pp[0]; pmax = pp[0];
    qmin = qq[0]; qmax = qq[0];
    for (i=1; i < ntxy; i++) {
	if (pp[i] < pmin) pmin=pp[i];
	if (qq[i] < qmin) qmin=qq[i];
	if (pp[i] > pmax) pmax=pp[i];
	if (qq[i] > qmax) qmax=qq[i];
    }
    dp = (pmax-pmin)/(npx-1);
    dq = (qmax-qmin)/(npy-1);

    wilson_init (ptx*10);

    for (iy=0; iy < npy; iy++) {
	q = (qmin + iy * dq) * dy / dt;
	ss->flt[5] = 0.5*q*(1.-q);
	ss->flt[6] = q*q-1.;
	ss->flt[7] = -0.5*q*(1.+q);
	
	for (ix=0; ix < npx; ix++) {
	    p = (pmin + ix * dp) * dx / dt;
	    f = p*p*(p*p-1.) + q*q*(q*q-1.);
	    s0 = 2. + 0.75*f;
	    ss->flt[0] = -f;
	    ss->flt[1] = 0.25*f;
	    ss->flt[2] = 0.5*p*(1.-p);
	    ss->flt[3] = p*p-1.;
	    ss->flt[4] = -0.5*p*(1.+p);

	    aa = sf_allocatehelix (ptx+1);
	    for (ia=0; ia < ptx+1; ia++) {
		aa->lag[ia] = ia+1;
		aa->flt[ia] = 0.;
	    }

	    a0 = wilson_factor (niter, 2.*s0, ss, aa, false, 1.e-6);
	    aa = compress (aa, eps);

	    for (ia=0; ia < aa->nh; ia++) {
		aa->flt[ia] = 0.;
	    }
 
	    sf_warning("%d %d %d", ix, iy, aa->nh);
	    a0 = wilson_factor (niter, 2.*s0, ss, aa, false, 1.e-6);
 
	    regrid(3, m, n, aa);

	    i = ix + npx * iy;
	    pfilt->hlx[i] = aa;
	}
    }

    for (i=0; i < ntxy; i++) {
	ix = 0.5 + (pp[i] - pmin)/ dp;
	iy = 0.5 + (qq[i] - qmin)/ dq;
	pfilt->pch[i] = ix + npx * iy;
    }

    sf_warning("Done with filters");

    sf_floatread (pp,ntxy,in);
    nhelicon_init (pfilt);
    nhelicon_lop (false,false,ntxy,ntxy,pp,qq);
    sf_floatwrite(qq,ntxy,out);

    exit(0);
}





