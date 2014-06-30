/* Bi-linear interpolation in 4-D for complex-valued data volume
NB: prepared for 5D interpolation using FFT
 */
/*
  Copyright (C) 2014 Xi'an Jiaotong University, UT Austin (Pengliang Yang)
  
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
#include <complex.h>
/*^*/

#include "lint4.h"

static int    nx,  ny, np, nq;
static float  ox,  oy, op, oq;
static float  dx,  dy, dp, dq;
static float *xx, *yy, *pp, *qq;

void lint4_init( int n1, float o1, float d1 /* first grid axis */, 
		 int n2, float o2, float d2 /* second grid axis */, 
		 int n3, float o3, float d3 /* third grid axis */,
		 int n4, float o4, float d4 /* fourth grid axis */,
		 float* x, float* y, float *p, float *q /* coordinates */) 
/*< initialize >*/
{
    nx = n1;
    ox = o1;
    dx = d1;
    xx = x;

    ny = n2;
    oy = o2;
    dy = d2;
    yy = y;

    np=n3;
    op=o3;
    dp=d3;
    pp=p;
    
    nq=n4;
    oq=o4;
    dq=d4;
    qq=q;
}

void lint4_lop( bool adj, bool add, int nm, int nd, sf_complex* mm, sf_complex *dd) 
/*< linear operator >*/
{
    int   ix, iy, ip, iq, im, id;
    float fx, fy, fp, fq;
    float gx, gy, gp, gq;
    float f;
    sf_complex h;

    if (nm != nx*ny*np*nq) sf_error("%s: wrong data size: %d != %d",__FILE__,nm,nx*ny*np*nq);
    sf_cadjnull (adj, add, nm, nd, mm, dd);
    
    for (id= 0; id< nd; id++) {

	f = (xx[id]-ox)/dx;
	ix = (int) f;  
	if ( ix < 0 || ix > nx-2) continue;
	fx=f-ix;
	gx= 1.-fx;
	
	f = (yy[id]-oy)/dy;
	iy = (int) f;  
	if ( iy < 0 || iy > ny-2) continue;
	fy=f-iy;
	gy= 1.-fy;
	
	f = (pp[id]-op)/dp;
	ip = (int) f;  
	if ( ip < 0 || ip > np-2) continue;
	fp=f-ip;
	gp= 1.-fp;
	
	f = (qq[id]-oq)/dq;
	iq = (int) f;  
	if ( iq < 0 || iq > nq-2) continue;
	fq=f-iq;
	gq= 1.-fq;
	
	im = ix+nx*iy+nx*ny*ip+nx*ny*np*iq;
	
	if (adj) {
	    mm[im  ]    += gx * gy * gp * gq * dd[id];
	    mm[im+1]    += fx * gy * gp * gq * dd[id];
	    mm[im+nx]   += gx * fy * gp * gq * dd[id];
	    mm[im+nx+1] += fx * fy * gp * gq * dd[id];

	    mm[im+nx*ny  ]    += gx * gy * fp * gq * dd[id];
	    mm[im+1+nx*ny]    += fx * gy * fp * gq * dd[id];
	    mm[im+nx+nx*ny]   += gx * fy * fp * gq * dd[id];
	    mm[im+nx+1+nx*ny] += fx * fy * fp * gq * dd[id];

	    mm[im+nx*ny*np  ]    += gx * gy * gp * fq * dd[id];
	    mm[im+1+nx*ny*np]    += fx * gy * gp * fq * dd[id];
	    mm[im+nx+nx*ny*np]   += gx * fy * gp * fq * dd[id];
	    mm[im+nx+1+nx*ny*np] += fx * fy * gp * fq * dd[id];

	    mm[im+nx*ny+nx*ny*np  ]    += gx * gy * fp * fq * dd[id];
	    mm[im+1+nx*ny+nx*ny*np]    += fx * gy * fp * fq * dd[id];
	    mm[im+nx+nx*ny+nx*ny*np]   += gx * fy * fp * fq * dd[id];
	    mm[im+nx+1+nx*ny+nx*ny*np] += fx * fy * fp * fq * dd[id];
	} else {
	  h=0.0;
	  h+=mm[im  ] * gx * gy * gp * gq;
	  h+=mm[im+1] * fx * gy * gp * gq;
	  h+=mm[im+nx]* gx * fy * gp * gq;
	  h+=mm[im+nx+1]* fx * fy * gp * gq;

	  h+=mm[im+nx*ny  ]  * gx * gy * fp * gq;
	  h+=mm[im+1+nx*ny]  * fx * gy * fp * gq;
	  h+=mm[im+nx+nx*ny]  * gx * fy * fp * gq;
	  h+=mm[im+nx+1+nx*ny] * fx * fy * fp * gq;

	  h+=mm[im+nx*ny*np  ] * gx * gy * gp * fq;
	  h+=mm[im+1+nx*ny*np] * fx * gy * gp * fq;
	  h+=mm[im+nx+nx*ny*np] * gx * fy * gp * fq;
	  h+=mm[im+nx+1+nx*ny*np]* fx * fy * gp * fq;

	  h+= mm[im+nx*ny+nx*ny*np  ] * gx * gy * fp * fq;
	  h+= mm[im+1+nx*ny+nx*ny*np] * fx * gy * fp * fq;
	  h+= mm[im+nx+nx*ny+nx*ny*np] * gx * fy * fp * fq;
	  h+= mm[im+nx+1+nx*ny+nx*ny*np]* fx * fy * fp * fq;

	  dd[id] += h;
	}
    }
}
