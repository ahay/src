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

#ifndef _norsar_h

struct point {
    float x;
    float z;
};
/*^*/

struct heptagon {
    struct point x0;
    float angle;
    float ampl;
    struct point x1;
    char cf;
    float srcAng;
};
/*^*/

struct grid {
    float *ampl;
    float *time;
    float *dirx;
    float *dirz;
    float *srcx;
    float *srcz;
    float *invgeo;
    int *flag;
};
/*^*/

#define IN_BOUNDS     	0
#define OUT_OF_BOUNDS 	1
#define CAUSTIC	      	2
#define END	      	3
#define NEWC	      	4
/*^*/

#define ROUND(x)	 ((int) (x))
/*^*/

#endif

#define DIST(A, B, C, D) sqrt((float) ((A)-(C))*((A)-(C))+((B)-(D))*((B)-(D)))
#define EPS 0.00001


/*----------------------------------------------------------------------------*/

float dist (struct point pt1, struct point pt2)
/*< distance >*/
{
    return DIST(pt1.x, pt1.z, pt2.x, pt2.z);
}

/*----------------------------------------------------------------------------*/

bool equal (float x, float X)
{
    return fabs(x - X) < EPS? true : false;
}

/*----------------------------------------------------------------------------*/

bool bigger (float x, float X)
/*< bigger >*/
{
    return (x - X) > EPS? true : false;
}


/*----------------------------------------------------------------------------*/

bool smaller (float x, float X)
/*< smaller >*/
{
    return (x - X) < -EPS? true : false;
}


/*----------------------------------------------------------------------------*/

struct point makepoint (float x, float z)
/*< make a point >*/
{
    struct point temp;

    temp.x = x;
    temp.z = z;
    return temp;
}

/*----------------------------------------------------------------------------*/

float slope (struct point pti, struct point ptj)
/*< Calculates the slope of the line that goes through points pti - ptj >*/
{
    return (ptj.z - pti.z) / (ptj.x - pti.x);
}

/*----------------------------------------------------------------------------*/

float ordinate (struct point pti, float m)
/*< Calculates the intersection of the line that goes through point
 * pti with a slope m with the z axis. >*/
{
    return pti.z - m * pti.x;
}


/*----------------------------------------------------------------------------*/

int bel (float X, float x1, float x2)
/*< Check to see whether x1 <= X <= x2 in case x1 < x2,
 * or x1 >= X >= x2 in case x1 >= x2. >*/
{
    return (x2>x1) ? ((X-x2<=EPS)&&(X-x1>=-EPS)) : ((X-x2>=-EPS)&&(X-x1<=EPS));
}

/*----------------------------------------------------------------------------*/

int eq_pts (struct point pt1, struct point pt2)
/*< equal points >*/
{
    float x, z;
    x = (float) fabs(pt1.x - pt2.x);
    z = (float) fabs(pt1.z - pt2.z);

    return ((x<=EPS)&&(z<=EPS)) ? 1 : 0;
}

/*----------------------------------------------------------------------------*/

void printcube (struct heptagon *cube, int nr, FILE *outfile)
/*< print cube for debugging >*/
{
    int ii;

    fprintf(outfile,"\nnr: %d", nr);
    for(ii=0;ii<nr;ii++) {
        fprintf(outfile,"\n\nx0: ( %f, %f)", cube[ii].x0.x, cube[ii].x0.z);
        fprintf(outfile,"\nangle: %f, flag: %d, amp: %f", cube[ii].angle * 180. /SF_PI, cube
[ii].cf, cube[ii].ampl);
        fprintf(outfile,"\nx1: ( %f, %f)", cube[ii].x1.x, cube[ii].x1.z);
    }
    return;
}
