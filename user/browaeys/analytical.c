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


struct pqvector
{
    float p[2];
    float q[2];
    float time;
    float sigma;
    int step;  /* step flag 0=undefined 1=dz 2=dx 3=dpz 4=dpx*/
};


#ifndef _analytical_h

typedef struct pqvector *hvec;
/*^*/

#endif


static int nx,nz;
static float dx,dz,ox,oz;
static float dpx,dpz;
static float **slow;
static float a[4];
static float b[4];


hvec hvec_init(float sigma /* evolution variable */,
               float time  /* traveltime */,
               float x     /* x position */,
               float z     /* depth position */,
               float px    /* px position */,
               float pz    /* pz position */)
/*< initialize phase space vector object >*/
{
    hvec pqvec;

    pqvec = (hvec)malloc(sizeof(struct pqvector));
    pqvec->step = 0;
    pqvec->sigma = sigma;
    pqvec->time = time;
    pqvec->q[0] = z;
    pqvec->q[1] = x;
    pqvec->p[0] = pz;
    pqvec->p[1] = px;

    return (pqvec);
}


void nc4_init(void)
/*< initialize Candi and Neri algorithm coefficients >*/
{
    a[0] = ( 2. + pow(2,1./3.) + pow(2,-1./3.) )/6.;
    a[1] = ( 1. - pow(2,1./3.) - pow(2,-1./3.) )/6.;
    a[2] = a[1];
    a[3] = a[0];

    b[0] = 0.;
    b[1] = 1./(2.-pow(2,1./3.));
    b[2] = 1./(1.-pow(2,2./3.));
    b[3] = b[1];
}


void slowg_lininterp(float *ssg, float x, float z)
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
	ssg[0] += (slow[izm+1][ixm]-slow[izm][ixm])/dz*(1.-tz);
	ssg[0] += (slow[izm+1][ixm+1]-slow[izm][ixm+1])/dz*tz;

        /* slowness x-gradient linear interpolation */
	ssg[1] += (slow[izm][ixm+1]-slow[izm][ixm])/dx*(1.-tx);
	ssg[1] += (slow[izm+1][ixm+1]-slow[izm+1][ixm])/dx*tx;

    }

    return;
}


float slow_bilininterp(float x, float z)
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

	ss += slow[izm][ixm]*(1.-tx)*(1.-tz);
	ss += slow[izm][ixm+1]*tx*(1.-tz);
	ss += slow[izm+1][ixm+1]*tx*tz;
	ss += slow[izm+1][ixm]*(1.-tx)*tz;       

	}

    return (ss);
}


hvec nc4_sigmastep(hvec pqvec, float ds)
/*< 4th order symplectic 2-D algorithm (Neri and Candy) marching in sigma >*/
{
    int i;
    float qi[2], pi[2], ssg[2], ss;

    for (i=0; i < 4; i++) {

        /* update for next step */
	pi[0] = pqvec->p[0]; /* pz */
	pi[1] = pqvec->p[1]; /* px */

	qi[0] = pqvec->q[0]; /* z */
	qi[1] = pqvec->q[1]; /* x */

        /* slowness interpolations in q space */
	ss = 0.;
	ss = slow_bilininterp(qi[1],qi[0]);

        /* slowness gradient interpolations in q space */
	ssg[0] = 0.;
	ssg[1] = 0.;
	slowg_lininterp(ssg,qi[1],qi[0]);

        /* advance p */
	pqvec->p[0] = pi[0] + b[i]*ss*ssg[0]*ds;
	pqvec->p[1] = pi[1] + b[i]*ss*ssg[1]*ds;

        /* advance q */
	pqvec->q[0] = qi[0] + a[i]*ds*(pqvec->p[0]);
	pqvec->q[1] = qi[1] + a[i]*ds*(pqvec->p[1]);

        /* advance traveltime (equivalent to Liouville transport equation) */
	pqvec->time += b[i]*ss*ss*ds;

    }
    pqvec->sigma += ds;

    return (pqvec);
}


float nc4_cellstep(hvec pqvec)
/*< sigma step from phase space cells step >*/
{
    float ds, dsz, dsx, dspz, dspx, ssg[2];

    /* Step size to exit from bottom or top */
    dsz = dz/fabs(pqvec->p[0]);

    /* Step size to exit from sides */
    dsx = dx/fabs(pqvec->p[1]);

    /* Step sizes to exit from slowness cell */ 
    ssg[0] = 0.;
    ssg[1] = 0.;

    slowg_lininterp(ssg,pqvec->q[1],pqvec->q[0]);

    dspz = dpz/fabs(ssg[0]);
    dspx = dpx/fabs(ssg[1]);

    /* Minimum sigma step size */

    /* step dz */
    ds = dsz;
    pqvec->step = 1;

    if (dsx < ds) {
        /* step dx */
	ds = dsx;
	pqvec->step = 2;
    }
    if (dspz < ds){
        /* step dpz */
	ds = dspz;
	pqvec->step = 3;
    }
    if (dspx < ds){
        /* step dpx */
	ds = dspx;
	pqvec->step = 4;
    }

    return (ds);

}

