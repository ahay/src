/* 2-D Laplacian operator */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

#include "laplacian.h"
#include "wilson.h"
#include "compress.h"

static int type,n1,n2,n12;
static float d1,d2, dd1, dd2, center, corner, *work1=NULL, *work2=NULL, **vt=NULL;
static sf_tris tri1=NULL,tri2=NULL;

void laplacian_init(int type1          /* operator type */,
		    int nz, int nx     /* dimensions */,
		    float dz, float dx /* sampling */,
		    float **vt1 /* [nx][nz] (v*t)^2 */)
/*< initialize >*/
{
    int i;
    float s1, s2, b0;
    const float tol=1.e-6;
    sf_filter ss, bb;

    type = type1;
    n1 = nz;
    n2 = nx;
    n12 = n1*n2;
    d1 = 1./(dz*dz);
    d2 = 1./(dx*dx);

    switch(type) {
	case 0:
	    center = -2.0*(d1+d2);
	    break;
	case 1:
	    corner = (d1+d2)/12.0;
	    s1 = (5*d1-d2)/6.0;
	    s2 = (5*d2-d1)/6.0;
	    d1 = s1;
	    d2 = s2;
	    center = -2.0*(2.0*corner+d1+d2);
	    break;
	case 2:
	    tri1 = sf_tridiagonal_init(n1);
	    sf_tridiagonal_const_define(tri1,10.0/(12.0*d1),1.0/(12.0*d1),false);
	    work1 = sf_floatalloc(n1);

	    tri2 = sf_tridiagonal_init(n2);
	    sf_tridiagonal_const_define(tri2,10.0/(12.0*d2),1.0/(12.0*d2),false);
	    work2 = sf_floatalloc(n2);
	    break;
	case 3:
	    corner = 3.0*(d1/d2+d2/d1);

	    ss = sf_allocatehelix(4);
	    ss->flt[0] = ss->flt[2] = (10.0 + corner)/144.0;
	    ss->flt[1] = ss->flt[3] = (2 - corner)/288.0;
	    ss->lag[0] = 1;
	    ss->lag[1] = n1-1;
	    ss->lag[2] = n1;
	    ss->lag[3] = n1+1;

	    bb = sf_allocatehelix(n1+1);
	    for (i=0; i <= n1; i++) {
		bb->lag[i] = i+1;
		bb->flt[i] = 0.0;
	    }

	    wilson_init(n1*10);
	    b0 = wilson_factor(100, (50.0 - corner)/72.0, ss, bb, true, tol);
	    wilson_close();
	    sf_deallocatehelix(ss);

	    bb = compress(bb,tol);
	    sf_warning("nb=%d",bb->nh);

	    work1 = sf_floatalloc(n12);

	    corner = (d1+d2)/(12.0*b0*b0);
	    s1 = (5*d1-d2)/(6.0*b0*b0);
	    s2 = (5*d2-d1)/(6.0*b0*b0);
	    d1 = s1;
	    d2 = s2;
	    center = -2.0*(2.0*corner+d1+d2);

	    sf_polydiv_init(n1*n2,bb);
	    break;
	case 4:
	    dd1 = -d1/12.0;
	    dd2 = -d2/12.0;
	    d1 *= 4.0/3.0;
	    d2 *= 4.0/3.0;
	    center = -2.0*(d1+d2+dd1+dd2);
	    break;
	case 5:
	    vt = vt1;
	    break;
	default:
	    sf_error("%s: Unknown Laplacian type",__FILE__);
    }
}

void laplacian_close(void)
/*< free allocated storage >*/
{
    switch(type) {
	case 2:
	    sf_tridiagonal_close(tri1);
	    free(work1);

	    sf_tridiagonal_close(tri2);
	    free(work2);
	    break;
	case 3:
	    sf_polydiv_close();
	    free(work1);
	    break;
	default:
	    break;
    }
}

void laplacian(float **uin  /* [nx][nz] */, 
	       float **uout /* [nx][nz] */)
/*< apply >*/
{
    int i1, i2;
    float s1, s2, v;

    for (i2=0; i2 < n2; i2++) {
	uout[i2][0]    = 0.;
	uout[i2][n1-1] = 0.;
    }
    for (i1=0; i1 < n1; i1++) {
	uout[0][i1]    = 0.;
	uout[n2-1][i1] = 0.;
    }
    
    switch(type) {
	case 0:
    	    for (i2=1; i2 < n2-1; i2++) {
		for (i1=1; i1 < n1-1; i1++) {
		    uout[i2][i1] = 
			d1*(uin[i2][i1-1]+uin[i2][i1+1]) +
			d2*(uin[i2-1][i1]+uin[i2+1][i1]) +
			center*uin[i2][i1];
		}
	    }
	    break;
	case 1:
	    for (i2=1; i2 < n2-1; i2++) {
		for (i1=1; i1 < n1-1; i1++) {
		    uout[i2][i1] = 
			d1*(uin[i2][i1-1]+uin[i2][i1+1]) +
			d2*(uin[i2-1][i1]+uin[i2+1][i1]) +
			corner*(uin[i2+1][i1-1]+uin[i2+1][i1+1]  +
				uin[i2-1][i1-1]+uin[i2-1][i1+1]) +
			center*uin[i2][i1];
		}
	    }
	    break;
	case 2:
	    for (i1=0; i1 < n1; i1++) {
		work2[0] = 0.;
		for (i2=1; i2 < n2-1; i2++) {
		    work2[i2] = uin[i2-1][i1]+uin[i2+1][i1]-2.0*uin[i2][i1];
		}
		work2[n2-1] = 0.;
		sf_tridiagonal_solve(tri2,work2);
		for (i2=0; i2 < n2; i2++) {
		    uout[i2][i1] = work2[i2];
		} 
	    }
	    for (i2=0; i2 < n2; i2++) {
		work1[0] = 0.;
		for (i1=1; i1 < n1-1; i1++) {
		    work1[i1] = uin[i2][i1-1]+uin[i2][i1+1]-2.0*uin[i2][i1];
		}
		work1[n1-1] = 0.;
		sf_tridiagonal_solve(tri1,work1);
		for (i1=0; i1 < n1; i1++) {
		    uout[i2][i1] += work1[i1];
		} 
	    }
	    break;
	case 3:
	    for (i2=1; i2 < n2-1; i2++) {
		for (i1=1; i1 < n1-1; i1++) {
		    uout[i2][i1] = 
			d1*(uin[i2][i1-1]+uin[i2][i1+1]) +
			d2*(uin[i2-1][i1]+uin[i2+1][i1]) +
			corner*(uin[i2+1][i1-1]+uin[i2+1][i1+1]  +
				uin[i2-1][i1-1]+uin[i2-1][i1+1]) +
			center*uin[i2][i1];
		}
	    }
	    sf_polydiv_lop(false,false,n12,n12,uout[0],work1);
	    sf_polydiv_lop(true, false,n12,n12,uout[0],work1);
	    break;
	case 4:
	    for (i2=2; i2 < n2-2; i2++) {
		for (i1=2; i1 < n1-2; i1++) {
		    uout[i2][i1] = 
			dd1*(uin[i2][i1-2]+uin[i2][i1+2]) +
			dd2*(uin[i2-2][i1]+uin[i2+2][i1]) +
			d1*(uin[i2][i1-1]+uin[i2][i1+1]) +
			d2*(uin[i2-1][i1]+uin[i2+1][i1]) +
			center*uin[i2][i1];
		}
	    }
	    break;
	case 5:
	    for (i2=2; i2 < n2-2; i2++) {
		for (i1=2; i1 < n1-2; i1++) {
		    v = vt[i2][i1];
		    dd1 = (v*d1-1.0)*d1/12.0;
		    dd2 = (v*d2-1.0)*d2/12.0;
		    s1 = (4.0-v*d1)*d1/3.0;
		    s2 = (4.0-v*d2)*d2/3.0;
		    center = -2.0*(s1+s2+dd1+dd2);
		    
		    uout[i2][i1] = 
			dd1*(uin[i2][i1-2]+uin[i2][i1+2]) +
			dd2*(uin[i2-2][i1]+uin[i2+2][i1]) +
			s1*(uin[i2][i1-1]+uin[i2][i1+1]) +
			s2*(uin[i2-1][i1]+uin[i2+1][i1]) +
			center*uin[i2][i1];
		}
	    }
	    break;
	default:
	    sf_error("%s: Unknown Laplacian type",__FILE__);
    }
}

