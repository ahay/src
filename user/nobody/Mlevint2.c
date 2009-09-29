/* Leveler inverse interpolation in 2-D. */
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

#include "helix.h"
#include "bound.h"
#include "createhelix.h"

#include "levint2.h"

int main(int argc, char* argv[])
{
    int nd, m[2], a[2], center[2], gap[2], nm, n2, id, nr, niter, warmup;
    float *rr, *dd, *xx, *yy, o1, d1, o2, d2, eps;
    float vmin, xmin, ymin, vmax, xmax, ymax, vrange, vaver;
    filter aa;
    bool cdstep;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    /* create model */
    if (!sf_getint ("m1",&m[0])) sf_error("Need m1=");
    if (!sf_getint ("m2",&m[1])) sf_error("Need m2=");
    /* number of bins */
    sf_putint(out,"n1",m[0]);
    sf_putint(out,"n2",m[1]);

    /* create filter */
    if (!sf_getint("a1",&a[0])) a[0]=5;
    if (!sf_getint("a2",&a[1])) a[1]=3;
    /* filter size */
    center[0] = a[0]/2;
    center[1] = 0;
    gap[0] = gap[1] = 0;
    aa = createhelix(2, m, center, gap, a);

    /* create model */
    nm = m[0]*m[1];
    nr = nm + aa->nh;
    rr = sf_floatalloc(nr);
    for (id=0; id < nr; id++) {
	rr[id] = 0.;
    }
    free(aa->flt);
    aa->flt = rr+nm;
    
    bound(2,m,m,a,aa);

    /* initialize filter with Laplacian */
    aa->flt[a[0]-1]   = -4.;
    aa->flt[a[0]-2]   = 1.;
    aa->flt[a[0]]     = 1.;
    aa->flt[2*a[0]-1] = 1.;
    
    /* initialize inversion */
    if (!sf_getint("niter",&niter)) niter=20;
    /* number of iterations */
    if (!sf_getint("warmup",&warmup)) warmup=20;
    /* number of initial iterations */
    if (!sf_getfloat("eps",&eps)) eps=1.;
    /* regularization parameter */
    if (!sf_getbool("cdstep",&cdstep)) cdstep=false;
    /* if y, use conjugate directions */

    if (!sf_histint(in,"n1",&nd)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2) || 3 != n2) sf_error("Need n2=3 in input");

    /* create filter */
    xx = sf_floatalloc(nd);
    yy = sf_floatalloc(nd);
    dd = sf_floatalloc(nd);

    sf_floatread (yy,nd,in);
    sf_floatread (xx,nd,in);
    sf_floatread (dd,nd,in);
    
    vmin = xmin = ymin = FLT_MAX;
    vmax = xmax = ymax = -FLT_MAX;
    for (id=0; id < nd; id++) {
	if (vmin > dd[id]) vmin=dd[id]; if (vmax < dd[id]) vmax=dd[id];
	if (xmin > xx[id]) xmin=xx[id]; if (xmax < xx[id]) xmax=xx[id];
	if (ymin > yy[id]) ymin=yy[id];	if (ymax < yy[id]) ymax=yy[id];
    }
    vaver  = (vmax+vmin)/2.;
    vrange = (vmax-vmin)/2.;
    for (id=0; id < nd; id++) {
	dd[id] = (vaver - dd[id])/vrange;
    }
    
    if (!sf_getfloat("o1",&o1)) o1=xmin;
    if (!sf_getfloat("o2",&o2)) o2=ymin;
    if (!sf_getfloat("d1",&d1)) d1=(xmax-xmin)/(m[0]-1);
    if (!sf_getfloat("d2",&d2)) d2=(ymax-ymin)/(m[1]-1);

    sf_putfloat(out,"o1",o1);
    sf_putfloat(out,"o2",o2);
    sf_putfloat(out,"d1",d1);
    sf_putfloat(out,"d2",d2);

    levint2 (cdstep, niter, warmup, nd, xx, yy, dd, 
	     o1, d1, m[0], o2, d2, m[1], aa, rr, eps);
    
    sf_floatwrite (rr,nm,out);

    exit(0);
}

/* 	$Id: Mlevint.c 1107 2005-04-11 16:52:19Z fomels $	 */
