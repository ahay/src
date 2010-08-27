/* 3-D velocity continuation using helix */
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

#include "velcon3.h"
#include "helify.h"

static float *r, *ld, *rd, *rhs, *w, *a, *b1, *b2, **pt, **flt;
static float v1, v2, dv, t0, dt, a0, cr, g2;
static const float g1 = 0.631974, b = 0.122996;
static int n, nt, nx, nv, inv;
static const int nf = 17, n1 = 8;
static sf_filter aa;

static void diffxx (const float *f, float *fxx)
/* helix laplacian */
{
    int i;

    for (i=0; i < nx+1; i++) {
	fxx[i] = 0.;
    }
    for (i=nx+1; i < n-nx-1; i++) {
	fxx[i] =  -cr * f[i] + 
	    g1 * (f[i   -1] + f[i   +1] + f[i-nx  ] + f[i+nx  ]) + 
	    g2 * (f[i-nx-1] + f[i+nx+1] + f[i-nx+1] + f[i+nx-1]); 
    }
    for (i=n-nx-1; i < n; i++) {
	fxx[i] = 0.;
    }
}

void velcon3_init (int inv1, float w1, float w2, float d0, 
		   int mt, int mx, int my, int mv, float dd, float dx)
/*< initialize >*/
{
    int it;

    n = mx*my; 
    nx=mx; 
    nt = mt; 
    nv = mv;

    v1 = w1; v2 = w2; inv = inv1;
    t0 = d0; dt = dd;
    
    r = sf_floatalloc(n);
    ld = sf_floatalloc(n);
    rd = sf_floatalloc(n);
    rhs = sf_floatalloc(n);
    pt = sf_floatalloc2(n,nt);

    dv = 0.5*(v1*v1-v2*v2)/nv; 
    a0=(4.*dx*dx)/(dt*dv);

    g2 = (1.-g1)*0.5;
    cr =  4.*(g1 + g2);

    flt = sf_floatalloc2(nf,nt);
    w = sf_floatalloc(nt);
    a = sf_floatalloc(nt);
    b1 = sf_floatalloc(nt);
    b2 = sf_floatalloc(nt);

    aa = sf_allocatehelix(nf);
    free (aa->flt);

    for (it=0; it < n1; it++) {
	aa->lag[it] = it+1;
    }
    for (it=n1; it < nf; it++) {
	aa->lag[it] = nx + 2 - nf + it;
    }
}

void velcon3_apply (bool adj, float **p1, float **p2)
/*< apply >*/
{
    int iv1, iv2, ivs, iv;
    int it1, it2, its, it, i;
    float t, b0, *l;

    if (adj) {
	iv1=0;  iv2=nv-1; ivs= 1; 
	it1=nt-1; it2= 1; its=-1;
	for (i=0; i < n; i++) {
	    for (it=0; it < nt; it++) {
		pt[it][i] = p2[i][it];
	    }
	}
    } else {
	iv1=nv-1; iv2= 0; ivs=-1; 
	it1= 1; it2=nt-1; its= 1; 
	for (i=0; i < n; i++) {
	    for (it=0; it < nt; it++) {
		pt[it][i] = p1[i][it];
	    }
	}	
    }

    for (it=it1; it != it2; it += its) { 
	t=t0+dt*it;

	switch(inv) {
	    case 0: /* Pseudounitary */
		t=sqrtf(t); a[it]=a0/t; break; 
	    case 1: /* Claerbout's */
		a[it]=a0;		break;   
	    default: /* True-amplitude */
		a[it]=a0/t; 	  t=1.;	break;   
	}
	b0 = b*a[it];
	if (b0 > t) b0 = 0.;
	b1[it]=b0+t;
	b2[it]=b0-t;
	aa->flt = flt[it];
	w[it] = helify (a[it],-b2[it],n1,nf,aa->flt);
    }

    for (iv=iv1; iv != iv2; iv += ivs) { 
	sf_warning("%d %d;",iv,iv2);

	for (i=0; i < n; i++) {
	    rhs[i] = 0.;
	}

	for (it=it1; it != it2; it += its) { 
	    aa->flt = flt[it]; 
	    l = pt[it]; 
	    diffxx(l,ld);

	    for (i=0; i < n; i++) {
		rhs[i] += a[it]*l[i] + b1[it]*ld[i];
	    }
	    sf_polydiv_init (n, aa);
	    sf_polydiv_lop (true, false, n, n, r, rhs);
	    sf_polydiv_lop (false,false, n, n, r, rhs);
	    for (i=0; i < n; i++) {
		r[i] = rhs[i]*w[it];
	    } 
	    diffxx(r,rd);
	    for (i=0; i < n; i++) {
		rhs[i] = a[it]*r[i] + b1[it]*rd[i] 
		    -    a[it]*l[i] - b2[it]*ld[i]; 
		l[i] = r[i];
	    }
	}
    }
     sf_warning(".");

    for (i=0; i < n; i++) {
	for (it=0; it < nt; it++) {
	    if (adj) {	
		p1[i][it] += pt[it][i];
	    } else {
		p2[i][it] += pt[it][i];
	    }
	}
    }
}
