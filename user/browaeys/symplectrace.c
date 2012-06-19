/* Neri and Candy 4th order symplectic algorithm and traveltime integration */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

#include "symplectrace.h"

#ifndef _symplectrace_h

typedef struct pqvector *pqv;
/* abstract data type */
/*^*/

#endif

struct pqvector {
    float p[2];
    float q[2];
    float time;
    int step;  /* step type is undefined(0), dz(1), dx(2), dp(3) */
};
/* concrete data type */

static float a[4];
static float b[4];

void nc4_init()
/*< initialize Neri and Candi algorithm coefficients >*/
{
    a[0] = ( 2. + pow(2,1./3.) + pow(2,-1./3.) )/6.;
    a[1] = ( 1. - pow(2,1./3.) - pow(2,-1./3.) )/6.;
    a[2] = a[1];
    a[3] = a[0];

    b[0] = 0.;
    b[1] = 1./(2.-pow(2,1./3.));
    b[2] = 1./(1.-pow(2,2./3.));
    b[3] = b[1];

    return;
}

void hvec_init(pqv pqvec,
               float time  /* traveltime */,
               float x     /* x position */,
               float z     /* depth position */,
               float px    /* px position */,
               float pz    /* pz position */)
/*< initialize phase space vector object >*/
{
    pqvec->step = 0;
    pqvec->time = time;
    pqvec->q[0] = z;
    pqvec->q[1] = x;
    pqvec->p[0] = pz;
    pqvec->p[1] = px;

    return;
}

void slowg_lininterp(float *ssg, float x, float z, float **slow, int nx, int nz, float dx, float dz, float ox, float oz)
/*< slowness gradient linear interpolation >*/
{
    int ixm, izm;
    float tz, tx;

    ssg[0] = 0.;
    ssg[1] = 0.;

    izm = floorf((z-oz)/dz);
    ixm = floorf((x-ox)/dx);

    if ( (ixm >= 0) && (izm >= 0) && (ixm < (nx-1)) && (izm < (nz-1)) ) {
	    
	tz = (z-izm*dz-oz)/dz;
	tx = (x-ixm*dx-ox)/dx;

        /* slowness z-gradient linear interpolation */
	ssg[0] += (slow[ixm][izm+1]-slow[ixm][izm])/dz*(1.-tz);
	ssg[0] += (slow[ixm+1][izm+1]-slow[ixm+1][izm])/dz*tz;

        /* slowness x-gradient linear interpolation */
	ssg[1] += (slow[ixm+1][izm]-slow[ixm][izm])/dx*(1.-tx);
	ssg[1] += (slow[ixm+1][izm+1]-slow[ixm][izm+1])/dx*tx;

    }

    return;
}

float slow_bilininterp(float x, float z, float **slow, int nx, int nz, float dx, float dz, float ox, float oz)
/*< slowness bilinear interpolation >*/
{
    int ixm, izm;
    float ss, tz, tx;

    ss = 0.;

    izm = floorf((z-oz)/dz);
    ixm = floorf((x-ox)/dx);

    if ( (ixm >= 0) && (izm >= 0) && (ixm < (nx-1)) && (izm < (nz-1)) ) {
	    
	tz = (z-izm*dz-oz)/dz;
	tx = (x-ixm*dx-ox)/dx;

	ss += slow[ixm][izm]*(1.-tx)*(1.-tz);
	ss += slow[ixm+1][izm]*tx*(1.-tz);
	ss += slow[ixm+1][izm+1]*tx*tz;
	ss += slow[ixm][izm+1]*(1.-tx)*tz;       

	}

    return (ss);
}

void nc4_sigmastep(pqv pqvec, float ds, float *ssi, float **slow, int nx, int nz, float dx, float dz, float ox, float oz)
/*< 4th order symplectic 2-D algorithm (Neri and Candy) marching in sigma >*/
{
    int i;
    float ssg[2], ss;

    for (i=0; i < 4; i++) {

	/* pqvec->p[0] is pz */
	/* pqvec->p[1] is px */
	/* pqvec->q[0] is z  */
	/* pqvec->q[1] is x  */

        /* slowness interpolations in q space */
	ss = slow_bilininterp(pqvec->q[1],pqvec->q[0],slow,nx,nz,dx,dz,ox,oz);

        /* slowness gradient interpolations in q space */
	slowg_lininterp(ssg, pqvec->q[1],pqvec->q[0],slow,nx,nz,dx,dz,ox,oz);

        /* slowness and gradient eno interpolation */
	/* sf_eno2_apply(cvel,kz,kx,0.,0.,&v1,g1,BOTH); */
	
        /* advance p */
	pqvec->p[0] += b[i]*ss*ssg[0]*ds;
	pqvec->p[1] += b[i]*ss*ssg[1]*ds;

        /* advance q */
	pqvec->q[0] += a[i]*ds*(pqvec->p[0]);
	pqvec->q[1] += a[i]*ds*(pqvec->p[1]);

        /* advance traveltime (i.e. Liouville transport equation) */
	pqvec->time += b[i]*ss*ss*ds;

    }

    *ssi = ss;

    return;
}

float nc4_cellstep(pqv pqvec, float **slow, int nx, int nz, float dx, float dz, float ox, float oz, float dpx, float dpz)
/*< signed sigma step from phase space cells step >*/
{
    float ds, dsz, dsx, dspx, ssg[2];

    dsz = SF_HUGE;
    dsx = SF_HUGE;

    /* dspz = SF_HUGE; */
    dspx = SF_HUGE;

    /* linear step size to exit from bottom or top */
    if (pqvec->p[0] != 0.0) dsz = dz/(pqvec->p[0]);

    /* linear step size to exit from sides */
    if (pqvec->p[1] != 0.0) dsx = dx/(pqvec->p[1]);

    /* linear step sizes to exit from slownesses cell */ 
    ssg[0] = 0.;
    ssg[1] = 0.;

    slowg_lininterp(ssg,pqvec->q[1],pqvec->q[0],slow,nx,nz,dx,dz,ox,oz);

    /* if (ssg[0] != 0.0) dspz = dpz/(ssg[0]); */
    if (ssg[1] != 0.0) dspx = dpx/(ssg[1]);

    /* select minimum sigma step size */

    /* step dz */
    ds = dsz;
    pqvec->step = 1;

    if (fabs(dsx) < fabs(ds)) {
        /* step dx */
	ds = dsx;
	pqvec->step = 2;
    }
    if (fabs(dspx) < fabs(ds)){
        /* step dpx */
	ds = dspx;
	pqvec->step = 3;
    }

    return (ds);

}

void value_exitlevel(pqv pqvec, float ss, int *step, float *x2, float *z2, float *p2, float *t)
/*< exiting values from computational step >*/
{
    *step = pqvec->step;
    *t = pqvec->time;

    *x2 = pqvec->q[1];
    *z2 = pqvec->q[0];
    *p2 = pqvec->p[1]/ss;

    return;
}
