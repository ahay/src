/* Interface for time-to-depth conversion */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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

#include "upgradt2d.h"
#include "t2diter.h"

static int nt, *f_0, *m_0, nt0, nx0;
static float dt0, dx0, ot0, ox0;
static float *t_0, *x_0, *v, *vd, *vdt, *vdx, *p_0;
static upgrad upg;
static float *temp1, *temp2;

float bilinear(const float *tbl, float t0, float x0);

void t2d_init(int dim  /* model dimension */,
	      int *n   /* model size */,
	      float *d /* model sampling */,
	      int n_t0, float d_t0, float o_t0,
	      int n_x0, float d_x0, float o_x0)
/*< initialize >*/
{
    int i;

    nt = 1;
    for (i=0; i < dim; i++) {
	nt = nt*n[i];
    }

    upg = upgrad_init(dim,n,d);

    temp1 = sf_floatalloc(nt);
    temp2 = sf_floatalloc(nt);

    nt0 = n_t0; dt0 = d_t0; ot0 = o_t0;
    nx0 = n_x0; dx0 = d_x0; ox0 = o_x0;
}

void t2d_caustic(float* x0   /* x0 */,
		 int* f0     /* f0 (modified in place) */, 
		 int* n, float* d /* model */,
		 float thres /* thresholding */)
/*< caustic region (2D) >*/
{
    int i, j;

    for (i=0; i < n[0]; i++) {
	for (j=0; j < n[1]; j++) {
	    if (j > 0) {
		if (x0[j*n[0]+i] <= x0[(j-1)*n[0]+i]) {
		    f0[j*n[0]+i] = 0;
		    continue;
		}
		if ((x0[j*n[0]+i]-x0[(j-1)*n[0]+i]) > thres*d[1]) {
		    f0[j*n[0]+i] = 0;
		    continue;
		}
	    }
	    if (j < n[1]-1) {
		if (x0[(j+1)*n[0]+i] <= x0[j*n[0]+i]) {
		    f0[j*n[0]+i] = 0;
		    continue;
		}
		if ((x0[(j+1)*n[0]+i]-x0[j*n[0]+i]) > thres*d[1]) {
		    f0[j*n[0]+i] = 0;
		    continue;
		}
	    }
	}
    }
}

void t2d_set(float* t0 /* t0 */,
	     float* x0 /* x0 */,
	     int* f0   /* f0 */,
	     float* v0   /* slowness squared */, 
	     float* vd0  /* Dix velocity */,
	     float* vdt0 /* Dix derivative in t0 */, 
	     float* vdx0 /* Dix derivative in x0 */,
	     int* m0   /* mask */,
	     float* p0 /* preconditioner */)
/*< set operator >*/
{
    t_0 = t0; x_0 = x0; f_0 = f0;
    v = v0; vd = vd0;
    vdt = vdt0; vdx = vdx0;
    m_0 = m0; p_0 = p0;

    upgrad_set(upg,t_0,x_0,f_0);
}

void t2d_close(void)
/*< free >*/
{
    upgrad_close(upg);
}

void t2d_cost(float* cost)
/*< evaluate cost >*/
{
    int i;

    upgrad_forw(upg,x_0,cost,false);

    for (i=0; i < nt; i++) {
	cost[i] -= bilinear(vd,t_0[i],x_0[i])*bilinear(vd,t_0[i],x_0[i])*v[i];
    }

    /* mask out points of one-sided stencil */
    for (i=0; i < nt; i++) {
	if (f_0[i] >= 0) cost[i] = 0.;
    }

    /* extra mask */
    if (m_0 != NULL) {
	for (i=0; i < nt; i++) {
	    if (m_0[i] >= 0) cost[i] = 0.;
	}
    }
}

void t2d_dt(const float* ds /* perturbation */,
	    float* dt       /* linear prediction of dt */)
/*< linear predict dt >*/
{
    int i;

    upgrad_solve(upg,ds,temp1,NULL,true);
    
    for (i=0; i < nt; i++) {
	dt[i] = 0.5*temp1[i];
    }
}

void t2d_dx(const float* dt /* dt */,
	    float* dx       /* linear prediction of dx */)
/*< linear predict dx >*/
{
    int i;

    upgrad_forw(upg,dt,temp2,false);

    for (i=0; i < nt; i++) {
	temp2[i] = -temp2[i];
    }

    upgrad_solve(upg,temp2,dx,NULL,true);
}

void t2d_prec(bool adj, bool add, int nd, int nr, float *d, float *r)
/*< preconditioner >*/
{
    int i;

    sf_adjnull(adj,add,nd,nr,d,r);

    if (adj) {
	for (i=0; i < nt; i++) {
	    d[i] = p_0[i]*r[i];
	}
    } else {
	for (i=0; i < nt; i++) {
	    r[i] = p_0[i]*d[i];
	}
    }
}

void t2d_oper(bool adj, bool add, int nx, int nr, float *x, float *r)
/*< linear operator >*/
{
    int i;

    sf_adjnull(adj,add,nx,nr,x,r);

    if (adj) { /* known misfit compute update */

	/* extra mask */
	if (m_0 != NULL) {
	    for (i=0; i < nt; i++) {
		if (m_0[i] >= 0) r[i] = 0.;
	    }
	}

	/* mask out points of one-sided stencil */
	for (i=0; i < nt; i++) {
	    if (f_0[i] >= 0) r[i] = 0.;
	}

	for (i=0; i < nt; i++) {
	    x[i] += bilinear(vd,t_0[i],x_0[i])*bilinear(vd,t_0[i],x_0[i])*r[i];
	}

	for (i=0; i < nt; i++) {
	    temp1[i] = bilinear(vd,t_0[i],x_0[i])*bilinear(vdt,t_0[i],x_0[i])*v[i]*r[i];
	}
	upgrad_inverse(upg,temp2,temp1,NULL,true);
	for (i=0; i < nt; i++) {
	    x[i] += temp2[i];
	}

	for (i=0; i < nt; i++) {
	    temp1[i] = -bilinear(vd,t_0[i],x_0[i])*bilinear(vdx,t_0[i],x_0[i])*v[i]*r[i];
	}
	upgrad_inverse(upg,temp2,temp1,NULL,true);
	upgrad_adj(upg,temp1,temp2,false);
	upgrad_inverse(upg,temp2,temp1,NULL,true);
	for (i=0; i < nt; i++) {
	    x[i] += temp2[i];
	}

	upgrad_adj(upg,temp1,r,false);
	upgrad_inverse(upg,temp2,temp1,NULL,true);
	upgrad_adj(upg,temp1,temp2,false);
	upgrad_inverse(upg,temp2,temp1,NULL,true);
	for (i=0; i < nt; i++) {
	    x[i] += temp2[i];
	}
    } else { /* known update compute misfit */

	for (i=0; i < nt; i++) {
	    r[i] += bilinear(vd,t_0[i],x_0[i])*bilinear(vd,t_0[i],x_0[i])*x[i];
	}

	upgrad_solve(upg,x,temp1,NULL,true);
	for (i=0; i < nt; i++) {
	    r[i] += bilinear(vd,t_0[i],x_0[i])*bilinear(vdt,t_0[i],x_0[i])*v[i]*temp1[i];
	}

	upgrad_forw(upg,temp1,temp2,false);
	upgrad_solve(upg,temp2,temp1,NULL,true);
	for (i=0; i < nt; i++) {
	    r[i] += -bilinear(vd,t_0[i],x_0[i])*bilinear(vdx,t_0[i],x_0[i])*v[i]*temp1[i];
	}

	upgrad_forw(upg,temp1,temp2,false);
	for (i=0; i < nt; i++) {
	    r[i] += temp2[i];
	}

	/* mask out points of one-sided stencil */
	for (i=0; i < nt; i++) {
	    if (f_0[i] >= 0) r[i] = 0.;
	}

	/* extra mask */
	if (m_0 != NULL) {
	    for (i=0; i < nt; i++) {
		if (m_0[i] >= 0) r[i] = 0.;
	    }
	}
    }
}

float bilinear(const float *tbl /* table input */,
	       float tin, float xin /* target position */)
/* bilinear interpolation */
{
    int neart, nearx;
    float distt, distx, val;

    /* find nearest grid */
    neart = (tin-ot0)/dt0;
    nearx = (xin-ox0)/dx0;

    if (neart < 0 || neart > nt0-1) sf_error("Interpolation beyond t0 scope.");
    if (nearx < 0 || nearx > nx0-1) sf_error("Interpolation beyond x0 scope.");

    if (neart == nt0-1) {
	neart = nt0-2; distt = 1.;
    } else {
	distt = (tin-(ot0+((float)neart)*dt0))/dt0;
    }
    
    if (nearx == nx0-1) {
	nearx = nx0-2; distx = 1.;
    } else {
	distx = (xin-(ox0+((float)nearx)*dx0))/dx0;
    }

    /* interpolation */
    val = (1.-distt)*(1.-distx)*tbl[nearx*nt0+neart]
	 +distt     *(1.-distx)*tbl[nearx*nt0+neart+1]
	 +(1.-distt)*distx     *tbl[(nearx+1)*nt0+neart]
	 +distt     *distx     *tbl[(nearx+1)*nt0+neart+1];
    
    return val;
}
