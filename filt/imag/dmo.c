#include <math.h>

#include <rsf.h>

#include "dmo.h" 
#include "halfint.h"

static bool inv;
static int nt2,nt,nx,ntx,mint,n,type;
static float vel,t0,dt,dx,dh,h,h2,d1,tmax, **tmp;

void dmo_init (float vel1, bool inv1, float t01, float dt1, float dx1,
	       int nt1, int nx1, float h1, int mint1, int n1, int type1)
{
    vel = vel1;
    inv = inv1;
    t0 = t01;
    dt = dt1;
    dx = dx1;
    nt = nt1;
    nx = nx1;
    h = h1;
    mint = mint1-1;
    n = n1;
    type = type1;

    nt2 = nt;
    if (nt%2) nt2++;

    ntx = nt*nx;
    h2=h*h;
    d1=fabsf(h)/n; 
    dh=d1*d1;
    tmax=1./(t0+(nt-1)*dt);

    tmp = sf_floatalloc2(nt2,nx);
}

void dmo_close(void)
{
    free (tmp[0]);
    free (tmp);
}

void dmo_lop (bool adj, bool add, int n1, int n2, float *dat1, float *dat2)
{
    float fx,ft,gx,gt,ama,x1,rat,x0,sl,sm,z,t,bx,amp, b1,r;
    int ix,iz,x,it,i1,i,kx,n0, id1, id2;

    sf_adjnull(adj,add,n1,n2,dat1,dat2);

    halfint_init(!adj,true,nt2,1.-1./nt2);
    
    if (!adj) {
	for (ix=0; ix < nx; ix++) {
	    for (it=0; it < nt; it++) { 
		tmp[ix][it] = dat1[it+ix*nt];
	    }
	    if (nt%2) tmp[ix][nt] = 0.;
	    halfint (tmp[ix]);
	}
    } else {
	for (i=0; i < nt2*nx; i++) { 
	    tmp[0][i] = 0.;
	}
    }

    for (iz=mint; iz < nt; iz++) { 
	z = t0 + dt * iz;
	x0 = 0.;
	n0 = 1;

	for (it=nt-1; it > iz; it--)  {    
	    t = t0 + dt * it; 
	    rat=(z*z)/(t*t); 
	    b1=h2*rat;

	    x0=sqrtf(h2-b1); 
	    n0=1+x0/d1;

	    sm=x0/b1;  
	    sl=t*sm*d1; 
	    if(sl < dt) break;

	    if (tmax > vel*sm) {
		r = z*tmax; 
	    } else {
		r = z*vel*sm;
	    }
	    if (r >= 1.) continue;

	    amp=(dt/sl)*sqrtf(z/(nt*dt));

	    switch (type) {
		case 0:
		    amp *= sqrtf(dh/b1);
		    if (adj && inv) amp *= h2*(h2+x0*x0)/(b1*b1);
		    break;
		case 1:
		    amp *= sqrtf(h2*dh)/b1;
		    if (adj && inv) amp *= (h2+x0*x0)/b1;
		    break;
		case 2:
		    amp *= (h2/b1)*sqrtf(dh/b1);
		    if (adj && inv) amp *= (h2+x0*x0)/h2;
		    break;
		case 3: default:
		    amp *= sqrtf(h2*(h2+x0*x0)/b1)/b1;
		    break;
	    }

	    bx=x0/dx; 
	    fx=fabsf(bx);

	    ix = fx; 
	    fx = fx-ix; 
	    gx = 1.-fx; 
	    x=ix+1;
	    if(x >= nx) continue;

	    for (kx=0; kx < nx-x; kx++) {
		id1 = it + (kx+ix)*nt;
		id2 = it + (kx+ x)*nt;

		if(adj) {
		    tmp[kx][iz] += amp*(dat2[id1]*gx + dat2[id2]*fx);
		} else {
		    ama=tmp[kx][iz]*amp;
		    dat2[id1] += gx*ama;
		    dat2[id2] += fx*ama;
		}
	    }

	    for (kx=x; kx < nx; kx++) {
		id1 = it + (kx-ix)*nt;
		id2 = it + (kx- x)*nt;

		if(adj) {
		    tmp[kx][iz] += amp*(dat2[id1]*gx + dat2[id2]*fx);
		} else {
		    ama=tmp[kx][iz]*amp;
		    dat2[id1] += gx*ama;
		    dat2[id2] += fx*ama;
		}
	    }
        }
    

	for (i1=0; i1 < n0; i1++) {  
	    x1= x0-d1*i1;       
	    b1=h2-x1*x1; 
	    sm=x1/b1;

	    if (tmax > vel*sm) {
		r = z*tmax;
	    } else {
		r = z*vel*sm;
	    }
	    if(r >= 1.) continue;

	    bx=x1/dx; 
	    fx=fabsf(bx);
	    ix = fx; 
	    fx = fx-ix; 
	    gx = 1.-fx; 
	    x=ix+1;        
	    if(x >= nx) continue;

	    t=z/sqrtf(b1/h2);
	    ft = (t-t0)/dt; 
	    it = ft; 
	    ft = ft-it; 
	    it = it+1; 
	    i=it+1; 
	    gt=1.-ft;
	    if(it >= nt) continue;

	    amp=sqrtf(z/(nt*dt));

	    switch (type) {
		case 0:
		    amp *= sqrtf(dh/b1);
		    if (adj && inv) amp *= h2*(h2+x1*x1)/(b1*b1);
		    break;
		case 1:
		    amp *= sqrtf(h2*dh)/b1;
		    if (adj && inv) amp *= (h2+x1*x1)/b1;
		    break;
		case 2:
		    amp *= (h2/b1)*sqrtf(dh/b1);
		    if (adj && inv) amp *= (h2+x1*x1)/h2;
		    break;
		case 3: default:
		    amp *= sqrtf(h2*(h2+x1*x1)/b1)/b1;
	    }

	    for (kx=0; kx < nx-x; kx++) {
		id1 = (kx+ix)*nt;
		id2 = (kx+ x)*nt;

		if(adj) {
		    tmp[kx][iz] += amp*(gt*gx*dat2[it+id1] +
					ft*gx*dat2[i +id1] +
					gt*fx*dat2[it+id2] +
					ft*fx*dat2[i +id2]);
		} else {
		    ama= tmp[kx][iz]*amp;
		    dat2[it+id1] += gt*gx*ama;
		    dat2[i +id1] += ft*gx*ama;
		    dat2[it+id2] += gt*fx*ama;
		    dat2[i +id2] += ft*fx*ama;
		}
	    }

	    for (kx=x; kx < nx; kx++) {
		id1 = (kx-ix)*nt;
		id2 = (kx- x)*nt;

		if(adj) {
		    tmp[kx][iz] += amp*(gt*gx*dat2[it+id1] +
					ft*gx*dat2[i +id1] +
					gt*fx*dat2[it+id2] +
					ft*fx*dat2[i +id2]);
		} else {
		    ama= tmp[kx][iz]*amp;
		    dat2[it+id1] += gt*gx*ama;
		    dat2[i +id1] += ft*gx*ama;
		    dat2[it+id2] += gt*fx*ama;
		    dat2[i +id2] += ft*fx*ama;
		}
	    }
	}
    }

    if (adj) {
        for (ix=0; ix < nx; ix++) {
	    if (nt%2) tmp[ix][nt] = 0.;
	    halfint (tmp[ix]);
	    
	    for (it=0; it < nt; it++) {
		dat1[it+ix*nt] += tmp[ix][it];
	    }
	}
    }

    halfint_close();
}
