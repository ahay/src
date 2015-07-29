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
#include "interp.h"

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

void mark_pts_outofbounds (struct heptagon *cube, int nr, float x0, float x1, float z0, float z1)
/*< "mark_pts_outofbounds" raises a flag on any point of the wavefront 
* that goes out of bounds. The boundary is define by points x0, x1,
* z0, z1. >*/
{
    int ii;
    struct point pt;
	
    for(ii=0;ii<nr;ii++) {
	pt = cube[ii].x1;
	if(smaller(pt.x, x0) || bigger(pt.x, x1) || smaller(pt.z, z0) || bigger(pt.z, z1)) 
	    cube[ii].cf=OUT_OF_BOUNDS;
    }
    return;
}

void makeup (struct heptagon *cube, int *nr)
/*< "makeup" takes off, from the wavefront, any point that goes
 * out of boundaries or that belongs to a caustic. >*/
{
    int ii;
    struct heptagon temp;
    int sp;

    sp = 0;
    for(ii=0;ii<*nr;ii++) {
	temp = cube[ii];
	switch (temp.cf) {
	    case OUT_OF_BOUNDS:
		if(cube[(ii-1+*nr)%*nr].cf != OUT_OF_BOUNDS) 		{
		    if(cube[(ii-sp-2+*nr)%*nr].cf!=END) 
			cube[(ii-sp-1+*nr)%*nr].cf = END;
		    else sp++;						}
	    case CAUSTIC:
		sp++;
		break;
	    case IN_BOUNDS:
	    default:
		cube[(ii-sp+*nr)%*nr] = temp;
		break;
	}
    }
    *nr -= sp;
    return;
}

void amplitudes (struct heptagon *cube,int  nr)
/*< This subroutine obtains the amplitudes of the new wavefront
* by calculating the geometrical spreading factor and 
* multiplying it with the old amplitude value.
*
* Notice that in order to calculate the amplitude for ending
* points of the wavefront; i.e., ii==0 or ii==nr-1, the code
* "wraps" to the other end of the wavefront by using moduli
* nr. >*/ 
{
    int ii;
    struct point A0, A1, A2;
    struct point B0, B1, B2;
    float R1, R2, r1, r2;

    for(ii=0;ii<nr;ii++) {
	A1 = cube[ii].x0;
	B1 = cube[ii].x1;

	if(cube[(ii-1+nr)%nr].cf==END) {
	    A0 = A1; B0 = B1;
	} else {
	    A0 = cube[(ii-1+nr)%nr].x0;
	    B0 = cube[(ii-1+nr)%nr].x1;
	}

	if(cube[(ii+nr)%nr].cf==END) {
	    A2 = A1; B2 = B1;
	} else {
	    A2 = cube[(ii+1+nr)%nr].x0;
	    B2 = cube[(ii+1+nr)%nr].x1;
	}

	r1 = dist(A0, A1); r2 = dist(A1, A2);
	R1 = dist(B0, B1); R2 = dist(B1, B2);
	if(R1+R2 > 0)
	    cube[ii].ampl *= sqrt((r1+r2) / (R1+R2));
    }
    return;
}

void interpolation (struct heptagon *cube, int *nr, int nrmax, float DSmax)
/*< interpolation checks the distance between two contiguos points in
* the wavefront, if its bigger than DSmax, then it calculates how
* many rays does it has to interpolate in between.
*
* interpolation calls a pair of subroutines which lay in ./interp.c
* Both subroutines use a third order polynome interpolation. >*/
{
    int ii, jj;
    struct point pt0, pt1, pt2, pt3;
    float A0, A1, A2, A3, s;
    float an0, an1, an2, an3;
    int sp, nnc;
    struct heptagon temp;

    for(ii=0;ii<*nr;ii++) {

	if(cube[ii].cf == END)  continue;	
	nnc = ROUND(dist(cube[ii].x1, cube[(ii+1+*nr)%*nr].x1) / DSmax);
	if(nnc==0) continue;

        if(*nr + 1 >= nrmax)
            sf_error("The wavefront has grown bigger than the amount of \nmemory you alocated for it. %d", *nr);

	pt1 = cube[ii].x1;
	pt2 = cube[(ii+1+*nr)%*nr].x1;
	A1 = cube[ii].ampl;
        A2 = cube[(ii+1+*nr)%*nr].ampl;
	an1 = cube[ii].angle;
	an2 = cube[(ii+1+*nr)%*nr].angle;
	if(ii==*nr-1) an2 += 2.*SF_PI;
	/* cf2 = cube[(ii+1+*nr)%*nr].cf; */

	sp=ii-1;
	while(cube[(sp+*nr)%*nr].cf==NEWC) sp--;
	sp = (sp+*nr)%*nr;
	if(cube[sp].cf==END) {
	    pt0 = pt1;
	    A0 = A1;
	    an0 = an1;
	} else {
	    pt0 = cube[sp].x1;
	    A0 = cube[sp].ampl;
	    an0 = cube[sp].angle;
	    if(ii==0) an0 -= 2.*SF_PI;
	}

	sp=ii+2;
	while(cube[(sp+*nr)%*nr].cf==NEWC) sp++;
	sp = (sp+*nr)%*nr;
	if(cube[sp].cf==END) {
	    pt3 = pt2;
	    A3 = A2;	
	    an3 = an2;
	} else {
	    pt3 = cube[sp].x1;
	    A3 = cube[sp].ampl;
	    an3 = cube[sp].angle;
	    if(ii>=*nr-2) an3 += 2.*SF_PI;
	}

	for(jj=0;jj<nnc;jj++) {	
	    s = (float) (jj + 1) / (nnc + 1);
	    temp.x0 = cube[ii].x0;
	    temp.cf = NEWC;
/* Call interpolating subroutines to obtain ray angle, amplitude, and
   position on the wavefront, for the new points.		*/
	    temp.angle = realinterp (pt0, pt1, pt2, pt3, an0, an1, an2, an3, s);
	    temp.ampl = realinterp (pt0, pt1, pt2, pt3, A0, A1, A2, A3, s);
	    temp.x1 = ptinterp (pt0, pt1, pt2, pt3, s);
/* Push new wavefront point into cube				*/
	    push(temp, cube, nr, ++ii);
	}
    }
    return;
}
