/* Traveltime and amplitude estimation using wavefront construction. */
/*
  Copyright (C) 1993 The Board of Trustees of Stanford University
  
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

#include "norgen.h"
#include "raytrace.h"

#define LN2      0.69314718

static int N, inter, nx, nz;
static float TETAMAX, alpha2, ox, oz, dx, dz;
static struct point length;

static int s(float x, float z);

#define CHECKX(X)       (X)<ox ? ox : (X)>ox+length.x ? ox+length.x :(X)
#define CHECKZ(Z)       (Z)<oz ? oz : (Z)>oz+length.z ? oz+length.z :(Z)

void raytrace_init(float tetamax, int n, float alpha, int inter1,
		   int nx1, int nz1, float ox1, float oz1, float dx1, float dz1,
		   struct point length1)
/*< initialize >*/
{
    TETAMAX=tetamax;
    N=n;
    alpha2=alpha;
    inter=inter1;
    nx=nx1; nz=nz1;
    ox=ox1; oz=oz1;
    dx=dx1; dz=dz1;
    length=length1;
}

/*----------------------------------------------------------------------------*/

void lomax (struct point pt0, float angle, float *vel, 
	    float dt, int nt, float T, float *temp)
/*< Lomax's method >*/
{
    int ii, nn;
    float x, z, x1, z1, x2, z2;
    float Xv, Zv, Xn, Zn;
    float Xv1, Zv1, Xn1, Zn1;
    float Xv2, Zv2, Xn2, Zn2;
    float V, V1, V2, cv, cn;
    float P, W, cw, sw;
    float a1, a2;

    a1 = TETAMAX * T / (4 * SF_PI * N);
    a2 = -4 * LN2 * TETAMAX * TETAMAX / (alpha2 * N * N);

    W = 1;
    for(ii=1; ii<=N; ii++) W += 2 * exp(a2 * ii * ii);

    x = x1 = x2 = pt0.x;
    z = z1 = z2 = pt0.z;

    for(ii=0;ii<nt;ii++) { 
	cw = cos(angle); sw = sin(angle);
	Xv = Xn = x;
	Zv = Zn = z;
	V = velo (Xv, Zv, vel);
	P = a1 * V;
	Xv1 = (Xn1 = Xv + P * cw);
	Zv1 = (Zn1 = Zv - P * sw);
	Xv2 = (Xn2 = Xv - P * cw);
	Zv2 = (Zn2 = Zv + P * sw);
	V1 = velo (Xv1, Zv1, vel);
	V2 = velo (Xv2, Zv2, vel);

	for(nn=1; nn<=N; nn++) {
	    Xn1 = Xn; Zn1 = Zn;
	    Xn = Xn2; Zn = Zn2;
	    Xv2 = Xv; Zv2 = Zv;
	    Xv = Xv1; Zv = Zv1;
	    cv = velo (Xv, Zv, vel);
	    cn = velo (Xn, Zn, vel);

	    Xv1 += a1 * cv * cw;
	    Zv1 -= a1 * cv * sw;
	    Xn2 -= a1 * cn * cw;
	    Zn2 += a1 * cn * sw;

	    V += exp(a2*nn*nn)*(cv + cn);
	    V1 += exp(a2*nn*nn)*(velo (Xv1, Zv1, vel) + velo (Xn1, Zn1, vel));
	    V2 += exp(a2*nn*nn)*(velo (Xv2, Zv2, vel) + velo (Xn2, Zn2, vel));

	}

	V /= W;
	V1 /= W;
	V2 /= W;
	x += V * dt * sw;
	z += V * dt * cw;
	angle += atan((V2 - V1) * dt / (2*P));
    }
    temp[0] = angle;
    temp[1] = x;
    temp[2] = z;

    return;
}

/*----------------------------------------------------------------------------*/

float velo (float X, float Z, float *vel)
/*< velocity >*/
{
    float u1, u2, xo, zo;
    float V0, V1, V2, V3;

    if(!inter)
	return vel[s(X, Z)];

    xo = (float) ROUND((X - ox) / dx) * dx + ox;
    zo = (float) ROUND((Z - oz) / dz) * dz + oz;

    u1 = (X - xo) / dx;
    u2 = (Z - zo) / dz;

    V0 = vel[s(xo             , zo)];
    V1 = vel[s(xo + dx, zo)];
    V2 = vel[s(xo             , dz + zo)];
    V3 = vel[s(xo + dx, dz + zo)];

    return V0*(1-u2)*(1-u1) + V1*u1*(1-u2) + V2*(1-u1)*u2 + V3*u1*u2;
}

/*----------------------------------------------------------------------------*/

static int s(float x, float z)
{
    x = CHECKX(x);
    z = CHECKZ(z);

    return ROUND((x-ox)/dx)*nz + ROUND((z-oz)/dz);
}
