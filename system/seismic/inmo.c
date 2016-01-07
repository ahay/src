/* Inverse normal moveout. */
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

#include <math.h>

#include <rsf.h>

#include "stretch4.h"
#include "bandpass.h"

static int nh, nt, ns, jump;
static float *stretch, **dense2, **sparse2; 
static map4* nmo;
static sf_bands spl;

void inmo_init(const float* vel /* velocity */,
	       const float* off /* offset */,
	       int nh1 /* number of offsets */,
	       float h0, float dh, int CDPtype, int ix,
	       int nt1 /* time samples */, bool slow,
	       float t0, float dt, float eps, bool half,
	       int jump1 /* subsampling */)
/*< initialization >*/
{
    int it, ih;
    float f, h, *coord;
    const int nw=4;

    nh = nh1;
    ns = nt1;
    jump = jump1;
    nt = (ns-1)*jump+1;
    dt = dt/jump;

    nmo = (map4*) sf_alloc(nh,sizeof(map4));
    stretch = sf_floatalloc(nt);

    for (ih=0; ih < nh; ih++) {
	nmo[ih] = stretch4_init (nt, t0, dt, nt, eps);
	
	h = off[ih] + (dh/CDPtype)*(ix%CDPtype); 
	if (half) h *= 2;
	h = h*h - h0*h0;
	    
	for (it=0; it < nt; it++) {
	    f = t0 + it*dt;
	    if (slow) {
		f = f*f + h*vel[it]*vel[it];
	    } else {
		f = f*f + h/(vel[it]*vel[it]);
	    }
	    if (f < 0.) {
		stretch[it]=t0-10.*dt;
	    } else {
		stretch[it] = sqrtf(f);
	    }
	}

	stretch4_define (nmo[ih],stretch);
    }

    coord = sf_floatalloc(nt);
    for (it=0; it < nt; it++) {
	coord[it] = it;
    }

    spl = sf_spline_init(nw,ns);
    sf_int1_init (coord, 0.0, jump, ns, sf_spline_int, nw, nt, 0.0);
    free(coord);
    
    dense2  = sf_floatalloc2(nt,nh);
    sparse2 = sf_floatalloc2(ns,nh);
}

void inmo_close(void)
/*< free allocated storage >*/
{
    int ih;

    for (ih=0; ih < nh; ih++) {
	stretch4_close(nmo[ih]);
    }
    free(nmo);
    free(stretch);
     
    sf_int1_close();
    sf_banded_close(spl);

    free(*dense2);
    free(dense2);
    free(*sparse2);
    free(sparse2);
}

void inmo(const float *trace   /* input trace [nt] */,
	  float **gather /* output CMP gather [nt][nh] */)
/*< apply inverse nmo >*/
{
    int ih;

    for (ih = 0; ih < nh; ih++) {
	stretch4_apply (false,nmo[ih],(float*)trace,gather[ih]);
    }
}


void nmostack(float **gather /* input CMP gather */,
	      float *trace /* output trace */)
/*< apply NMO and stack >*/
{
    int ih, it;

    for (it=0; it < nt; it++) {
	trace[it] = 0.0f;
    }

    for (ih = 0; ih < nh; ih++) {
	stretch4_invert (false,nmo[ih],stretch,gather[ih]);
	for (it=0; it < nt; it++) {
	    trace[it] += stretch[it];
	}
    }

    /* normalize */
    for (it=0; it < nt; it++) {
	trace[it] /= nh;
    }
}


void subsample(float **dense, float **sparse) 
/*< subsample a gather >*/
{
    int ih, it, is;

    for (ih = 0; ih < nh; ih++) {
	for (it=is=0; it < nt; it+=jump, is++) {
	    sparse[ih][is] = dense[ih][it];
	}
    }
}


void interpolate(float **sparse, float **dense)
/*< interpolate to a denser sampling >*/
{
    int ih;

    for (ih = 0; ih < nh; ih++) {
	sf_banded_solve(spl,sparse[ih]);
	sf_int1_lop (false,false,ns,nt,sparse[ih],dense[ih]);
    }

}


void inmo_oper(int nx, const float *trace, float *trace2, void* user_data)
/*< operator to invert by GMRES >*/
{
    int it;
    
    /* forward operator */
    inmo(trace,dense2);
    subsample(dense2,sparse2);
    /* backward operator */
    interpolate(sparse2,dense2);
    nmostack(dense2,trace2);

    /* I + S (BF - I) */
    for (it=0; it < nt; it++) {
	trace2[it] -= trace[it];
    }
    bandpass(trace2);
    for (it=0; it < nt; it++) {
	trace2[it] += trace[it];
    }
}
	      
