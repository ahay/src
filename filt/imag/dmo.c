#include <math.h>

#include <rsf.h>

#include "dmo.h" 
#include "halfint.h"

static bool inv;
static int nt,nx,ntx,mint,n,type;
static float vel,t0,dt,dx,dh,h,h2,d1,tmax, *tmp;

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

    ntx = nt*nx;
    h2=h*h;
    d1=fabsf(h)/n; 
    dh=d1*d1;
    tmax=1./(t0+(nt-1)*dt);

    tmp = sf_floatalloc(ntx);
}

void dmo_close(void)
{
    free (tmp);
}

void dmo_lop (bool adj, bool add, int n1, int n2, dat1, dat2)
{
    float fx,ft,gx,gt,ama,x1,rat,x0,sl,sm,z,t,bx,amp, b1,r;
    int ix,iz,maxt,x,it,i1,i,kx,n0, im, id1, id2;

    sf_adjnull(adj,add,n1,n2,dat1,dat2);

    halfint_init(!adj,true,nt,1.-1./nt);
    
    if (!adj) {
	for (i=0; i < ntx; i++) { 
	    tmp[i] = dat1[i];
	}
        for (ix=0; ix < nx; ix++) {
	    halfint (tmp+nt*ix);
	}
    } else {
	for (i=0; i < ntx; i++) { 
	    tmp[i] = 0.;
	}
    }

    for (iz=mint; iz < nt; iz++) { 
	z = t0 + dt * iz;

	for (it=nt-1; it > iz, it--)  {    
	    t = t0 + dt * it; 
	    rat=(z*z)/(t*t); 
	    b1=h2*rat;

	    x0=sqrtf(h2-b1); 
	    n0=1+x0/d1;

	    sm=x0/b1;  
	    sl=t*sm*d1; 
	    if(sl < dt) break;

	    r = z*amax1(tmax,vel*sm);               
	    if (r >= 1.) continue;

	    amp=(dt/sl)*sqrt(z/(nt*dt));

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
	    if(x >= nx) cycle;

	    for (kx=0; kx < nx-x; kx++) {
		im = iz + kx*nt;
		id1 = it + (kx+ix)*nt;
		id2 = it + (kx+ x)*nt;

		if(adj) {
		    tmp[im] += amp*(dat2[id1]*gx + dat2[id2]*fx);
		} else {
		    ama=tmp[im]*amp;
		    dat2[id1] += gx*ama;
		    dat2[id2] += fx*ama;
		}
	    }

	    for (kx=x; kx < nx; kx++) {
		im = iz + kx*nt;
		id1 = it + (kx-ix)*nt;
		id2 = it + (kx- x)*nt;

		if(adj) {
		    tmp[im] += amp*(dat2[id1]*gx + dat2[id2]*fx);
		} else {
		    ama=tmp[im]*amp;
		    dat2[id1] += gx*ama;
		    dat2[id2] += fx*ama;
		}
	    }
        }
    

	for (i1=0; i1 < n0; i1++) {  
	    x1= x0-d1*i1;       
	    b1=h2-x1*x1; 
	    sm=x1/b1;

	    r = z*amax1(tmax,vel*sm);               
	    if(r >= 1.) cycle;

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
		im = iz + kx*nt;
		id1 = (kx+ix)*nt;
		id2 = (kx+ x)*nt;

		if(adj) {
		    tmp[im] += amp*(gt*gx*dat2[it+id1] +
				    ft*gx*dat2[i +id1] +
				    gt*fx*dat2[it+id2] +
				    ft*fx*dat2[i +id2]);
		} else {
		    ama= tmp[im]*amp;
		    dat2[it+id1] += gt*gx*ama;
		    dat2[i +id1] += ft*gx*ama;
		    dat2[it+id2] += gt*fx*ama;
		    dat2[i +id2] += ft*fx*ama;
		}
	    }

	    for (kx=x; kx < nx; kx++) {
		im = iz + kx*nt;
		id1 = (kx-ix)*nt;
		id2 = (kx- x)*nt;

		if(adj) {
		    tmp[im] += amp*(gt*gx*dat2[it+id1] +
				    ft*gx*dat2[i +id1] +
				    gt*fx*dat2[it+id2] +
				    ft*fx*dat2[i +id2]);
		} else {
		    ama= tmp[im]*amp;
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
	    halfint (tmp+ix*nt);
	    
	    for (i=0; i < ntx; i++) {
		dat1[i] += tmp[i];
	    }
	}
    }

    halfint_close();
}
