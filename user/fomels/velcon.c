/* Velocity continuation with finite differences */
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

static float *r, *ld, *rd, *rhs, **pt;
static float v1, v2, t0, dt, dv, a0;
static const float b=0.14867678; /* "one sixth" factor */
static int nt, nx, nv, inv;
static sf_tris slv;

static void diffxx (int nx, const float *f, float *fxx);

void velcon_init (int inv1 /* inversion type */, 
		  float v11, float v21 /* starting and ending velocities */, 
		  float t01            /* time origin */, 
		  int nt1              /* time samples */,
		  int nx1              /* midpoint samples */,
		  int nv1              /* velocity steps */,
                  float dt1            /* time sampling */,
		  float dx             /* midpoint sampling */)
/*< initialize >*/
{
    inv=inv1; v1=v11; v2=v21;
    t0=t01; nt=nt1; nx=nx1;
    nv=nv1; dt=dt1; 

    r = sf_floatalloc(nx);
    ld = sf_floatalloc(nx);
    rd = sf_floatalloc(nx);
    rhs = sf_floatalloc(nx);

    pt = sf_floatalloc2(nx,nt);

    dv = (v1-v2)/nv; 
    a0=(4.*dx*dx)/(dt*dv);
    slv = sf_tridiagonal_init(nx);
}

void velcon_close(void) 
/*< free allocated storage >*/
{
    free(r);
    free(ld);
    free(rd);
    free(rhs);
    free(*pt);
    free(pt);
    sf_tridiagonal_close(slv);
}

void velcon_lop (bool adj, bool add, int n1, int n2, float *p1, float *p2)
/*< linear operator >*/
{
    int iv1, iv2, ivs, iv;
    int it1, it2, its, it, ix;
    float v, t, a, b1, b2, diag, *l;

    sf_adjnull(adj,add,n1,n2,p1,p2);

    if (adj) {
	iv1=0;    iv2=nv; ivs= 1; 
	it1=nt-1; it2= 0; its=-1; 
    } else {
	iv1=nv-1; iv2=-1; ivs=-1;
	it1=1;    it2=nt; its= 1; 
    }
       
    for (it=0; it < nt; it++) {
	for (ix=0; ix < nx; ix++) {
	    pt[it][ix] = adj? p2[ix*nt+it]: p1[ix*nt+it];
	}
    }

    for (iv=iv1; iv != iv2; iv += ivs) {  
	v=sqrtf(v2+(iv+0.5)*dv); 
	for (ix=0; ix < nx; ix++) {
	    rhs[ix] = 0.;
	}
	for (it=it1; it != it2; it +=its) {  
	    t=t0+dt*it;
	    switch (inv) {
		case 0: /* Pseudounitary */
		    t=sqrtf(t)*v; 	
		    a=a0/t;
		    break;
		case 1: /* Claerbout's */ 
		    t *= v;	        
		    a=a0/v;
		    break;
		default: /* True-amplitude */
		    a=a0/(t*v); 	
		    t=v;
		    break;
	    }
	    b1=b*a+t; 
	    b2=b*a-t; 
	    diag=a-2*b2;
	    sf_tridiagonal_const_define(slv,diag,b2,true);

	    l = pt[it]; 		         
	    diffxx(nx,l,ld);
	    for (ix=0; ix < nx; ix++) {
		r[ix] = rhs[ix] + a*l[ix] + b1*ld[ix];
	    }
	    sf_tridiagonal_solve(slv,r);
	    diffxx(nx,r,rd);
	    for (ix=0; ix < nx; ix++) {
		rhs[ix] = a*r[ix] + b1*rd[ix] - a*l[ix] - b2*ld[ix];
		l[ix] = r[ix];
	    }
	}
    }

    for (it=0; it < nt; it++) {
	for (ix=0; ix < nx; ix++) {
	    if (adj) {
		p1[ix*nt+it] += pt[it][ix];
	    } else {
		p2[ix*nt+it] += pt[it][ix];
	    }
	}
    }
}

static void diffxx (int nx, const float *f, float *fxx) 
/* compute second finite difference */
{
    int ix;

    fxx[0] = f[1]-f[0];
    for (ix=1; ix < nx-1; ix++) {
	fxx[ix] = f[ix-1]-2*f[ix]+f[ix+1];
    }
    fxx[nx-1] = f[nx-2]-f[nx-1];
}


