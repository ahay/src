#include <math.h>

#include <rsf.h>

#include  "kirchnew.h"

static int nt, nx, sw;
static float t0, dt, dx, *vrms;

void kirchnew_init (float *vrms_in, float t0_in, float dt_in, float dx_in, 
		    int nt_in, int nx_in, int sw_in)
{
    vrms = vrms_in;
    t0 = t0_in;
    dt = dt_in;
    dx = dx_in;
    nt = nt_in;
    nx = nx_in;
    sw = sw_in;
}

void kirchnew_lop (bool adj, bool add, int nm, int nd, 
		   float *modl, float *data)
{
    int im, id1, id2, ix,iz,it,ib,iy, minx[2],maxx[2], is,i;
    float amp,t,z,b,db,f,g;

    sf_adjnull(adj,add,nm,nd,modl,data);

    maxx[0] = nx;	
    minx[1] = 1;
    for (iz=0; iz < nt-1; iz++) {     
	z = t0 + dt * iz;		/* vertical traveltime */
	for (it=nt-1; it >= iz+1; it--) { 
	    t = t0 + dt * it;		/* time shift */
	    b = sqrtf(t*t - z*z); 
	    db = dx*b*2./(vrms[iz]*t);
	    if(db < dt || sw == 1) break;

	    f = 0.5*vrms[iz]*b/dx; 
	    iy = f; f = f-iy; 
	    i = iy+1; g = 1.-f;

	    if(i >= nx)	continue;
	    amp = (z / (t+dt)) * sqrtf(nt*dt / (t+dt)) * (dt / db);

	    minx[0] = i;  
	    maxx[1] = nx-i;
	    for (is=0; is < 2; is++) {  
		iy = -iy; 
		i = -i;		/* two branches of hyperbola */
		for (ix=minx[is]; ix < maxx[is]; ix++) {
		    im = ix*nt+iz;
		    id1 = (ix+iy)*nt+it;
		    id2 = (ix+i)*nt+it;
		    if( adj) {	
			modl[im] += amp*(data[id1]*g + data[id2]*f);
		    } else {
			data[id1] += modl[im]*amp*g;
			data[id2] += modl[im]*amp*f;
		    }
		}
	    }
	}

	for (ib=0; ib < nx; ib++) {	   
	    b = dx*ib*2./vrms[iz]; /* space shift */ 
	    iy = ib; 
	    t = hypotf(z,b); 
	    db = dx*b*2./(vrms[iz]*t);
	    if(db > dt || sw == 2) break;

	    f = (t-t0)/dt; 
	    it = f; f = f-it; 
	    i = it+1; g = 1.-f;
	    if(it >= nt) break;

	    amp = (z / (t+dt)) * sqrtf(nt*dt / (t+dt)); 
	    if(ib == 0) amp *= 0.5;

	    minx[0] = iy; 
	    maxx[1] = nx-iy;
	    for (is=0; is < 2; is++) {	
		iy = -iy; /* two branches of hyperbola */
		for (ix=minx[is]; ix <	maxx[is]; ix++) {
		    im = ix*nt+iz;
		    id1 = (ix+iy)*nt+it;
		    id2 = (ix+iy)*nt+i;

		    if( adj) {	
			modl[im] += amp*(data[id1]*g + data[id2]*f);
		    } else {
			data[id1] += modl[im]*amp*g;
			data[id2] += modl[im]*amp*f;
		    }
		}
	    }
	}
    }
}
