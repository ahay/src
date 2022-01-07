/* Helper functions for Gauss-Seidel */
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

/* Phase space grid */
static int nz,nx,na;
static float oz,ox,oa;
static float dz,dx,da;
static float dzi,dxi,dai;
static const char code[] = "xatzl";

void gsray_init(int nz1, int nx1, int na1,
		float oz1, float ox1, float oa1,
		float dz1, float dx1, float da1)
/*< initialize >*/
{
    nz = nz1; oz = oz1; dz = dz1;
    nx = nx1; ox = ox1; dx = dx1;
    na = na1; oa = oa1; da = da1;

    dzi = -1.0/dz;
    dxi = -1.0/dx;
    dai = -1.0/da;
}

float gs_update(float ***t, 
		float Hp,
		float Hq,
		float da_dsigma,
		int iz, int ix, int ia, 
		/* float ss, float ssz, float ssx,  slowness and slowness gradient */
		/* float cs, float sn, cosine and sine of the phase angle 
		   float ss2, */
		float dt_dsigma,
		float dl_dsigma,
		int iq /* what to compute */)
/*< Gauss-Seidel update >*/
{
    float ts, fz, fx, fa, dd;

    float fzmax, fzmin; 
    float fxmax, fxmin; 
    float famax, famin;

    /* Downwind scheme coefficients  */
    fz = dzi*Hp;  /* cs*ss; */
    fx = dxi*Hq;  /* sn*ss; */
    /* fa = dai*(cs*ssx - sn*ssz);*/
    fa = dai*da_dsigma;
    
    fzmin = SF_MAX(-fz,0.);
    fzmax = SF_MAX(fz,0.);
    
    fxmin = SF_MAX(-fx,0.);
    fxmax = SF_MAX(fx,0.);
    
    famin = SF_MAX(-fa,0.);
    famax = SF_MAX(fa,0.);

    /* Diagonal term */
    dd = (fxmax + fxmin + fzmax + fzmin + famax + famin);

    ts = 0.0;

    /* z-direction */
    if (iz == 0) {

	ts += fzmax*t[ia][ix][1];

    } else if (iz == (nz-1)) {

	ts += fzmin*t[ia][ix][nz-2];

    } else {

	ts += fzmax*t[ia][ix][iz+1] + fzmin*t[ia][ix][iz-1];

    }

    /* x-direction */
    if (ix == 0) {

	ts += fxmax*t[ia][1][iz];

    } else if (ix == (nx-1)) {

	ts += fxmin*t[ia][nx-2][iz];

    } else {

	ts += fxmax*t[ia][ix+1][iz] + fxmin*t[ia][ix-1][iz];

    }

    /* escape angle is periodic on ]-PI,+PI] */
    if (ia == 0) {

	ts += famax*t[1][ix][iz] + famin*t[na-1][ix][iz];

    } else if (ia == (na-1)) {

	ts += famax*t[0][ix][iz] + famin*t[na-2][ix][iz];

    } else {

	ts += famax*t[ia+1][ix][iz] + famin*t[ia-1][ix][iz];

    }

    switch(code[iq]) {
	case 't':
	    ts += dt_dsigma; /* ss2; ss*ss */
	    break;
	case 'l':
	    ts += dl_dsigma; /* ss in anisotropic wave p. */
	    break;
    }

    return (ts/dd);
}

char gs_color(char ***c, 
	      int iz, int ix, int ia, 
	      float ss, float ssz, float ssx, /* slowness and slowness gradient */
	      float cs, float sn /* cosine and sine of the phase angle */)
/*< Gauss-Seidel "color" >*/
{
    char cz;

    float fz;

    float fzmax, fzmin; 

    /* Downwind scheme coefficients  */
    fz = dzi*cs*ss;
    
    fzmin = SF_MAX(-fz,0.);
    fzmax = SF_MAX(fz,0.);
    
    /* z-direction */
    if (fzmax > 0. && iz < nz-1) {
	cz = c[ia][ix][iz+1];
    } else if (fzmin > 0. && iz > 0) {
	cz = c[ia][ix][iz-1];
    } else {
	cz = '0';
    }

    return cz;
}


void boundary_mat(float ***t, int iz, int ix, int ia, int iq)
/*< Initialize outgoing rays for boundary conditions >*/
{
    switch(code[iq]) {
	case 'x':
	    t[ia][ix][iz] = ox + ix*dx;
	    break;
	case 'a':
	    /* convert to degrees */
	    t[ia][ix][iz] = (oa + ia*da)*180./SF_PI;
	    break;
	case 't':
	case 'l':
	    t[ia][ix][iz] = 0.0;
	    break;
	case 'z':
	    t[ia][ix][iz] = oz + iz*dz;
	    break;
    }

    return;
}
	 
