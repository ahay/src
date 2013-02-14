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
#include <stdio.h>

#include <rsf.h>

#include "norsar.h"

#include "norgen.h"
/*^*/

#include "raytrace.h"

#define IN_BOUNDS     	0

static int gnx, gnz;

void norsar_init(int gnx1, int gnz1,
		 float TETAMAX, int N, float alpha2, int inter, 
		 int nx, int nz, float ox, float oz, float dx, float dz,
		 struct point length)
/*< initialize >*/
{
    gnx = gnx1;
    gnz = gnz1;
    raytrace_init(TETAMAX,N,alpha2,inter,
		  nx,nz,ox,oz,dx,dz,length);
}

/*----------------------------------------------------------------------------*/

void initial (struct point pos, struct heptagon *cube, float *vel, 
	      float dt, int nt, float T, int lomx, int nr, struct grid *out)
/*<
* This subroutine initializes "cube", a structure array that contains:
* the location of the points over the wavefront at a time t,
* take-off angles from these points, the ending points over the
* wavefront at a time t + nt*dt and a flag (which tells if the ray
* belongs to a caustic or if its out of bounds).
*
* cube[].x0		position of source on wavefront at time t
*			cube[].x0 and cube[].x1 are point structures;
*			i.e., they an x and z components.
*			
* cube[].angle		arriving angle from cube[].x1.x, cube[].x1.z
* cube[].ampl		amplitude
* cube[].x1		position of arriving ray on wavefront at time t + dt
* cube[].cf		flag
*
* cube is rewritten every time a wavefront is propagated.
*
* The number of elements of cube needs to be, as big as the number 
* of points on the wavefront. Notice that there is no way to estimate
* this number a priori, so a rough estimate for a higher bound should
* be given at the beginning of the program (nrmax).
>*/
{
    int ii;
    float temp, phi; /* ra[3]; */

    temp = 10.; /*/sqrt((float) nang);*/

    for(ii=0;ii<nr;ii++) {
        phi = ii * 2. * SF_PI / nr;
        cube[ii].x0 =    pos;
        cube[ii].angle = phi;
	cube[ii].cf = IN_BOUNDS;
    }

    wavefront(cube, nr, vel, dt, 1, T, lomx);

    for(ii=0;ii<nr;ii++) {
        cube[ii].x0 =   cube[ii].x1;
        cube[ii].ampl = temp;
    }

/* Make use of this subroutine to empty contents of output array 
   For an explanation on of this output array, look in the file
   ./gridding							*/
    for(ii=0;ii<gnx*gnz;ii++) {
	out->ampl[ii] = 0.;
	out->time[ii] = 0.;
	out->flag[ii] = 0;
    }

    return;
}

/*----------------------------------------------------------------------------*/

void wavefront (struct heptagon *cube, int nr, float *vel, 
		float dt, int nt, float T, int lomx)
/*<
* This subroutine propagates the wavefront, i.e; computes points over 
* the next wavefront. wavefront uses subroutine lomax which prototype 
* is explained here:
*
*
* void lomax (point, take-off angle, velocity model, time step,
*		number of time steps, period, array that stores
*		values passed back by lomax)
*
* The subroutine "lomax" propagates a waveray of a given period,
* through the velocity model starting from a source location x,z
* at a given take-off angle. It computes the ending location x',z' 
* and its landing angle, after a nt*dt time. The ouput is given
* by the array ra[].
*
* ra[0] = arriving angle of ray at x', z'
* ra[1] = x' position of ray
* ra[2] = z' position of ray
*
>*/
{
    int ii;
    float phi;
    struct point pt0;
    float ra[3];

    /* FOR EVERY POINT ON THE WAVEFRONT, OBTAIN ITS LOCATION ON
    *  THE NEW WAVEFRONT AND ARRIVING ANGLE */

    for(ii=0;ii<nr;ii++) {
	phi = cube[ii].angle;
	pt0 = cube[ii].x0;

	if(lomx) lomax (pt0, phi, vel, dt, nt, T, ra);
	else sf_error("Use lomx=1, no other option is available");
	cube[ii].x1 = makepoint(ra[1], ra[2]);
	cube[ii].angle = ra[0];
    }
    return;
}
