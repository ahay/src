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
#include <stdio.h>
#include <stdlib.h>

struct pqvector
{
    float p[2];
    float q[2];
    float time;
    float sigma;
};

typedef struct pqvector *hvec;

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


hvec nc4_step(hvec pqvec, float ds, float **slow,
              const int nx, const int nz, const float ox, const float oz, const float dx, const float dz)
/*< 4th order symplectic 2-D algorithm (Neri and Candy) marching in sigma >*/
{
    int i, ixm, izm;
    float qi[2], pi[2], ssg[2], ss, tx, tz;

    for (i=0; i < 4; i++) {

        /* update for next step */
	pi[0] = pqvec->p[0]; /* pz */
	pi[1] = pqvec->p[1]; /* px */

	qi[0] = pqvec->q[0]; /* z */
	qi[1] = pqvec->q[1]; /* x */

        /* slowness and gradient interpolations in q space */
	ss = 0.;
	ssg[0] = 0.;
	ssg[1] = 0.;

	izm = floorf((qi[0]-oz)/dz);
	ixm = floorf((qi[1]-ox)/dx);

	if ( (ixm >= 0) && (izm >= 0) && (ixm < (nx-1)) && (izm < (nz-1)) ) {
	    
	    tz = (qi[0]-izm*dz-oz)/dz;
	    tx = (qi[1]-ixm*dx-ox)/dx;

            /* slowness bilinear interpolation */
	    ss += slow[izm][ixm]*(1.-tx)*(1.-tz);
	    ss += slow[izm][ixm+1]*tx*(1.-tz);
	    ss += slow[izm+1][ixm+1]*tx*tz;
	    ss += slow[izm+1][ixm]*(1.-tx)*tz;

            /* slowness z-gradient linear interpolation */
	    ssg[0] += (slow[izm+1][ixm]-slow[izm][ixm])/dz*(1.-tz);
	    ssg[0] += (slow[izm+1][ixm+1]-slow[izm][ixm+1])/dz*tz;

            /* slowness x-gradient linear interpolation */
	    ssg[1] += (slow[izm][ixm+1]-slow[izm][ixm])/dx*(1.-tx);
	    ssg[1] += (slow[izm+1][ixm+1]-slow[izm+1][ixm])/dx*tx;

	}

        /* advance p */
	pqvec->p[0] = pi[0] + b[i]*ss*ssg[0]*ds;
	pqvec->p[1] = pi[1] + b[i]*ss*ssg[1]*ds;

        /* advance q */
	pqvec->q[0] = qi[0] + a[i]*ds*(pqvec->p[0]);
	pqvec->q[1] = qi[1] + a[i]*ds*(pqvec->p[1]);

        /* advance traveltime */
	pqvec->time += b[i]*ss*ss*ds;

    }
    pqvec->sigma += ds;

    return (pqvec);

}

float dsigmaz(float dz, float pz)
/*< metrics transformation dz -> dsigma >*/
{
    return (dz/pz);
}

float dsigmax(float dx, float px)
/*< metrics transformation dx -> dsigma >*/
{
    return (dx/px);
}

float dsigmap(float dpx, float ssx)
/*< metrics transformation dpx -> dsigma >*/
{
    return (dpx/ssx);
}

float analytical(float x1, float z1, float v1, const float *g1, float *p1,
		 float x2, float z2, float v2, const float *g2, float *p2)
/*< return traveltime >*/
{
    float d[2], a, g[2], v, z, dist, grad, disc, time, num, den, r, r1;
    float gg1, gg2, aa;

    d[0] = z2-z1;
    d[1] = x2-x1;

    a = 1.;

    num = v2-g2[0]*d[0]-g2[1]*d[1] - v1;
    den = v1+g1[0]*d[0]+g1[1]*d[1] - v2;

    if (fabsf(den) < fabsf(num)) {
	r = den/num;
	r1 = 1.+r;
	if (fabsf(r1) > 0.) a = 1./r1;
    } else if (fabsf(num) < fabsf(den)) {
	r = num/den;
	r1 = 1.+r;
	if (fabsf(r1) > 0.) a = r/r1;
    }
    
    if (a < 1. && a >= 0.5) {
	a = 1.;
    } else if (a < 0.5 && a > 0.) {
	a = 0.;
    }

    g[0] = g2[0]+a*(g1[0]-g2[0]);
    g[1] = g2[1]+a*(g1[1]-g2[1]);

    v = v1+v2;
    dist = d[0]*d[0]+d[1]*d[1];
    grad = g[0]*g[0]+g[1]*g[1];

    disc = v*v-dist*grad;
    if (disc > 0.) {
	z = 4*dist/(v+sqrtf(disc));
    } else {
	z = 4*dist/v;
    }

    p1[0] = d[0]-0.25*z*g[0];
    p1[1] = d[1]-0.25*z*g[1];
    disc = 1./hypotf(p1[0],p1[1]);
    p1[0] *= disc;
    p1[1] *= disc;

    p2[0] = d[0]+0.25*z*g[0];
    p2[1] = d[1]+0.25*z*g[1];    
    disc = 1./hypotf(p2[0],p2[1]);
    p2[0] *= disc;
    p2[1] *= disc;
    
    if (z > 0.) {
	time = dist + z*z*grad/48.;

	if (a > 1. && a < 0.) {
	    gg1 = (g1[0]-g2[0])*g1[0]+(g1[1]-g2[1])*g1[1];
	    gg2 = (g1[0]-g2[0])*g2[0]+(g1[1]-g2[1])*g2[1];

	    aa = (a-1.)*a;

	    time += aa*z*z*
		(-2.5*gg2 - 
		 (a*((7  + 2*a*(19 + 5*a*(8*a*(10 + 3*(a-3)*a)-37)))*gg1 + 
		     (13 - 2*a*(64 + 5*a*(8*a*(10 + 3*(a-3)*a)-43)))*gg2))/2. 
		 + 30*aa*aa*
		 (a*(3 - 6*a + 4*a*a)*(gg1 - gg2) + gg2)*logf(a/(a-1)))/120.;
	}
	
	time /= sqrtf(z);
    } else {
	time = 0.;
    }

    return time;
}
