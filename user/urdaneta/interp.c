/* Interpolation. */
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

//TODO - testing
#include <stdio.h>

static float ptmult (struct point A1, struct point A2);
static struct point ptdiff(struct point A1, struct point A2);

float realinterp (struct point pt0, struct point pt1, struct point pt2, struct point pt3, 
		  float f0, float f1, float f2, float f3, float var)
/*< real interpolation >*/
{
    float a, b, c, d;
    float A0, A1, dA0, dA1;

    A0 = f1;
    A1 = f2;
    dA0 = (f2 - f0) / dist(pt2, pt0);
    dA1 = (f3 - f1) / dist(pt3, pt1);

    d = A0;
    c = dA0;
    b = 3. * (A1 - A0) - 2. * dA0 - dA1;
    a = A1 - b - c - d;

    return ((a * var + b) * var + c) * var + d;
}

struct point ptinterp (struct point pt0, struct point pt1, struct point pt2, struct point pt3, 
		       float var)
/*< point interpolation >*/
{
    struct point a, b, c, d; 
    struct point X0, X1, dX0, dX1;
    struct point ei, ek;
    struct point Temp;
    float temp;

    ei = makepoint (1., 0.);
    ek = makepoint (0., 1.);

    X0 = pt1;
    X1 = pt2;
    temp = dist(pt2, pt0);
    dX0.x = ptmult(ptdiff(pt2, pt0), ei) / temp;
    dX0.z = ptmult(ptdiff(pt2, pt0), ek) / temp;
    temp = dist(pt3, pt1);
    dX1.x = ptmult(ptdiff(pt3, pt1), ei) / temp;
    dX1.z = ptmult(ptdiff(pt3, pt1), ek) / temp;
 
    d = X0;
    c = dX0;
    b.x = 3. * (X1.x - X0.x) - 2. * dX0.x - dX1.x;
    b.z = 3. * (X1.z - X0.z) - 2. * dX0.z - dX1.z;
    a.x = X1.x - b.x - c.x - d.x;
    a.z = X1.z - b.z - c.z - d.z;

    Temp.x = d.x + var * (c.x + var * (b.x + var * a.x));
    Temp.z = d.z + var * (c.z + var * (b.z + var * a.z));

    return Temp;
}

/*----------------------------------------------------------------------------*/

void push (struct heptagon newhpt, struct heptagon *cube, int *nr, int var)
/*< push >*/
{
    int ii;

    (*nr)++;
    for(ii=*nr;ii>var;ii--) {
	cube[ii] = cube[ii-1];
    }
    cube[var] = newhpt;

    return;
}

/*----------------------------------------------------------------------------*/

static float ptmult (struct point A1, struct point A2)
{
    return A1.x * A2.x + A1.z * A2.z;
}

/*----------------------------------------------------------------------------*/

static struct point ptdiff(struct point A1, struct point A2)
{
    struct point temp;

    temp.x = A1.x - A2.x;
    temp.z = A1.z - A2.z;
    return temp;
}

/*----------------------------------------------------------------------------*/

void TwoD_interp (struct grid *out, int gnx, int gnz)
/*< TwoD_interp is a simple subroutine that only averages amplitude and trave-
* time values for receivers that have missing values. Later on these 
* subroutine could be improve using bicubic splines or filters. >*/
{
    int ii;
    float ampl, time,dirx,dirz;

    for(ii=0;ii<gnx*gnz;ii++) {
        if(out->flag[ii]==0) {



            if(ii%gnz==0) {

		time  = out->time[ii+1-gnz] + out->time[ii+1+gnz];
		time += out->time[ii+gnz] + out->time[ii+1] + out->time[ii-gnz];
		time /= 5.;
		ampl  = out->ampl[ii+1-gnz] + out->ampl[ii+1+gnz];
		ampl += out->ampl[ii+gnz] + out->ampl[ii+1] + out->ampl[ii-gnz];
		ampl /= 5.;

		dirx  = out->dirx[ii+1-gnz] + out->dirx[ii+1+gnz];
		dirx += out->dirx[ii+gnz] + out->dirx[ii+1] + out->dirx[ii-gnz];
		dirx /= 5.;
		dirz  = out->dirz[ii+1-gnz] + out->dirz[ii+1+gnz];
		dirz += out->dirz[ii+gnz] + out->dirz[ii+1] + out->dirz[ii-gnz];
		dirz /= 5.;
            } else if(ROUND(ii/gnz)==0) {

		time  = out->time[ii-1+gnz] + out->time[ii+1+gnz];
		time += out->time[ii+gnz] + out->time[ii+1] + out->time[ii-1];
		time /= 5.;
		ampl  = out->ampl[ii-1+gnz]+ out->ampl[ii+1+gnz];
		ampl += out->ampl[ii+gnz] + out->ampl[ii+1] + out->ampl[ii-1];
		ampl /= 5.;
		dirx  = out->dirx[ii-1+gnz]+ out->dirx[ii+1+gnz];
		dirx += out->dirx[ii+gnz] + out->dirx[ii+1] + out->dirx[ii-1];
		dirx /= 5.;
		dirz  = out->dirz[ii-1+gnz]+ out->dirz[ii+1+gnz];
		dirz += out->dirz[ii+gnz] + out->dirz[ii+1] + out->dirz[ii-1];
		dirz /= 5.;
            } else if(ii%gnz==gnz-1) {
	
		time  = out->time[ii-1-gnz] + out->time[ii-1+gnz];
		time += out->time[ii-gnz] + out->time[ii-1] + out->time[ii+gnz];
		time /= 5.;
		ampl  = out->ampl[ii-1-gnz] + out->ampl[ii-1+gnz];
		ampl += out->ampl[ii-gnz] + out->ampl[ii-1] + out->ampl[ii+gnz];
		ampl /= 5.;
		dirx  = out->dirx[ii-1-gnz] + out->dirx[ii-1+gnz];
		dirx += out->dirx[ii-gnz] + out->dirx[ii-1] + out->dirx[ii+gnz];
		dirx /= 5.;
		dirz  = out->dirz[ii-1-gnz] + out->dirz[ii-1+gnz];
		dirz += out->dirz[ii-gnz] + out->dirz[ii-1] + out->dirz[ii+gnz];
		dirz /= 5.;
            } else if(ROUND(ii/gnz)==gnx-1) {

		time  = out->time[ii-1-gnz] + out->time[ii+1-gnz];
		time += out->time[ii-gnz] + out->time[ii-1] + out->time[ii+1];
		time /= 5.;
		ampl  = out->ampl[ii-1-gnz] + out->ampl[ii+1-gnz];
		ampl += out->ampl[ii-gnz] + out->ampl[ii-1] + out->ampl[ii+1];
		ampl /= 5.;
		dirx  = out->dirx[ii-1-gnz] + out->dirx[ii+1-gnz];
		dirx += out->dirx[ii-gnz] + out->dirx[ii-1] + out->dirx[ii+1];
		dirx /= 5.;
		dirz  = out->dirz[ii-1-gnz] + out->dirz[ii+1-gnz];
		dirz += out->dirz[ii-gnz] + out->dirz[ii-1] + out->dirz[ii+1];
		dirz /= 5.;
            } else {
		if(out->flag[ii-1]!=0) {
			time = out->time[ii-1];
		    ampl = out->ampl[ii-1];
		    dirx = out->dirx[ii-1];
		    dirz = out->dirz[ii-1];
		} else if (out->flag[ii+1]!=0) {
			time = out->time[ii+1];
		    ampl = out->ampl[ii+1];
		    dirx = out->dirx[ii+1];
		    dirz = out->dirz[ii+1];
		} else if (out->flag[ii-gnz]!=0) {
			time = out->time[ii-gnz];
		    ampl = out->ampl[ii-gnz];
		    dirx = out->dirx[ii-gnz];
		    dirz = out->dirz[ii-gnz];
		} else if (out->flag[ii+gnz]!=0) {
			time = out->time[ii+gnz];
		    ampl = out->ampl[ii+gnz];
		    dirx = out->dirx[ii+gnz];
		    dirz = out->dirz[ii+gnz];
		} else {
		   time = 0.;
		   ampl = 0.;
		   dirx = 0.;
		   dirz = 0.;
		}
            }

	    out->time[ii] = time;
	    out->ampl[ii] = ampl;
	    out->dirx[ii] = dirx;
	    out->dirz[ii] = dirz;
        }
    }
    return;
}
/*----------------------------------------------------------------------------*/

void movwavf (struct heptagon *cube, int nr)
/*< move wavefront >*/
{
    int ii;

    for(ii=0;ii<nr;ii++) {
	cube[ii].x0 = cube[ii].x1;
	if(cube[ii].cf==NEWC) cube[ii].cf = IN_BOUNDS;
    }
    return;
}

/******************************** END *****************************************/
