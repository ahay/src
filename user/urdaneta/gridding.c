/* Gridding. */
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

#include "norgen.h"
#include "interp.h"
#include "raytrace.h"

#define NUMmaxCELLperQUAD       50
#define NUM_MAX_REC		50

static int cells_of_quad (struct point *pt, int *cll);
static int common_receivers (int *cll, int num, int *rec);
static int seg_cross_chck (struct point pt0, struct point pt1, struct point pp0, struct point pp1, 
			   struct point *pt);
static void falling_rec (int *rec, struct grid *out, int *num_rec, 
			 struct point pt0, struct point pt1, struct point pt2, struct point pt3);
static void interp_rec(int *rec, struct grid *out, int num_rec, 
		       struct heptagon c0, struct heptagon c1, struct heptagon c2, struct heptagon c3, 
		       float tim, float *vel, int first);
static int cell (float x, float z);

static int gnx, gnz;
static float gdx, gdz, gox, goz;
static FILE * outfile;

struct point sourcePosition;

void gridding_init(int gnx1, int gnz1,
		   float gdx1, float gdz1, float gox1, float goz1,
		   FILE * outfile1, struct point  sourceOrigin)
/*< initialize >*/
{
    gnx = gnx1;
    gnz = gnz1;
    gdx = gdx1; gdz = gdz1;
    gox = gox1; goz = goz1;
    outfile = outfile1;
    sourcePosition =sourceOrigin;
}

void gridding (struct heptagon *cube, int nr, struct grid *out, float DSmax, float tim, float *vel, int first)
/*< Interpolates travel times and amplitude values from the wavefronts to
* a regular grid.
* The ouput goes to the array out. out is a pointer to structure of
* pointers which has the form:
*
* out->time	Contains travel-time table
* out->ampl	Contains amplitude table
* out->dirx X cosine of the ray direction
* out->dirz Z cosine of the ray direction
* out->ampl	Contains amplitude table
* out->srcx X cosine of the ray direction at source
* out->srcz Z cosine of the ray direction at source
* out->invgeo is the inverse geometrical spreading (i.e R in homogenous medium)
*
* out->flag	Tells which receivers have or have not time and
*		amplitude values. >*/
{
    int ii;
    int cll[NUMmaxCELLperQUAD], num_cells, num_rec;
    int rec[NUM_MAX_REC];
    struct point pt[4], pt_crss;
    struct heptagon c0, c1, c2, c3;

    for(ii=0;ii<nr;ii++) {

	if(cube[ii].cf==END) continue;

	c1 = cube[ii];	c2 = cube[(ii+1+nr)%nr];
	
        /* Compute points of ray cell pt0-pt1-pt2-pt3 */
        pt[0] = c1.x0; pt[1] = c2.x0;
        pt[2] = c1.x1; pt[3] = c2.x1;

        /* Find to which cells does quadrangle pt0-pt1-pt3-pt2 belongs to */
        num_cells = cells_of_quad (pt, cll);

        /* Find receivers (grid points) that are inside the cells */
	num_rec = common_receivers (cll, num_cells, rec);
        /* Check if ray pt0 - pt2 crosses ray pt1 - pt3 */
        if(seg_cross_chck (pt[0], pt[2], pt[1], pt[3], &pt_crss)) {
	    continue;
/*
            Check if rec[] falls inside crossed ray cell 
	       which is divided into two triangles		
	    falling_rec (rec, out, &num_rec, pt[0], pt[1], pt_crss, pt_crss);
	    falling_rec (rec, out, &num_rec, pt_crss, pt_crss, pt[2], pt[3]);
*/
        } else {

            /* In case they don't intersect, get which receivers
            fall inside the ray cell                             */
            falling_rec (rec, out, &num_rec, pt[0], pt[1], pt[2], pt[3]);

        }

	/* If there are no receivers inside the ray cell go to next
	   ray cell else interpolate ray cell values to receivers */

	if(num_rec!=0) {
	    if(cube[(ii-1+nr)%nr].cf==END) {
		c0 = c1;
	    } else {
		c0 = cube[(ii-1+nr)%nr];
	    }
	    if(c2.cf==END) {
		c3 = c2;
	    } else {
		c3 = cube[(ii+2+nr)%nr];
	    }

	    interp_rec(rec, out, num_rec, c0, c1, c2, c3, tim, vel, first);
	}
    }
    return;
}

/*----------------------------------------------------------------------------*/

static int cells_of_quad (struct point *pt, int *cll)
{
    int ii, jj;
    int numcols, numrows;
    float xmax, zmax, xmin, zmin;

    xmin = xmax = pt[0].x;
    zmin = zmax = pt[0].z;

    for(ii=1;ii<4;ii++) {
	if (pt[ii].x > xmax) xmax = pt[ii].x;
	else if (pt[ii].x < xmin) xmin = pt[ii].x;
	if (pt[ii].z > zmax) zmax = pt[ii].z;
        else if (pt[ii].z < zmin) zmin = pt[ii].z;
    }

    numcols = ROUND(cell(xmax, zmin)/gnz) - ROUND(cell(xmin, zmin)/gnz) + 1;
    numrows = cell(xmin, zmax) - cell(xmin, zmin) + 1;

    if(numcols * numrows > NUMmaxCELLperQUAD)	{
	fprintf(stderr, "\np0.x %f, p0.z: %f", pt[0].x, pt[0].z);
	fprintf(stderr, "\np1.x %f, p1.z: %f", pt[1].x, pt[1].z);
	fprintf(stderr, "\np2.x %f, p2.z: %f", pt[2].x, pt[2].z);
	fprintf(stderr, "\np3.x %f, p3.z: %f", pt[3].x, pt[3].z);
	fprintf(stderr, "\nnumcols: %d, numrows: %d" ,numcols, numrows);
	sf_error("%s: Grid size too small.",__FILE__); 
    }

    cll[0] = cell(xmin, zmin);
    for(ii=0;ii<numcols;ii++) {
	for(jj=0;jj<numrows;jj++)
	    cll[jj+ii*numrows] = cll[0] + jj + ii*gnz;
    }
    return numrows * numcols;
}

/*----------------------------------------------------------------------------*/

static int common_receivers (int *cll, int num, int *rec)
{
    int ii;

    for(ii=0;ii<num;ii++)
	rec[ii] = cll[ii];
    return num;
}

/*----------------------------------------------------------------------------*/

static int cell (float x, float z)
{
    x = (x<gox) ? gox : (x>=gox+gdx*(gnx-1)) ? gox+gdx*(gnx-1) : x;
    z = (z<goz) ? goz : (z>=goz+gdz*(gnz-1)) ? goz+gdz*(gnz-1) : z;
    return ROUND((z-goz)/gdz) + ROUND((x-gox)/gdx)*gnz;
}

/*----------------------------------------------------------------------------*/

static void falling_rec (int *rec, struct grid *out, int *num_rec, 
			 struct point pt0, struct point pt1, struct point pt2, struct point pt3)
/* Check which points in array rec[] fall inside ray cell 
* pt0-pt1-pt3-pt2-pt0
*
* In case a point falls outside the ray cell, it is erase
* from the array. */
{
    int ii;
    int sp;
    struct point receiv, pt_useless;

    sp = 0;
    for(ii=0;ii<*num_rec;ii++) {
	receiv = makepoint(ROUND(rec[ii]/gnz)*gdx+gox,(rec[ii]%gnz)*gdz+goz);

	if(!eq_pts (pt2, pt3)) 
	    if(seg_cross_chck(receiv, pt0, pt2, pt3, &pt_useless)!=0) continue;
	if(seg_cross_chck(receiv, pt0, pt3, pt1, &pt_useless)!=0) continue;

	if(!eq_pts (pt2, pt3))
	    if(seg_cross_chck(receiv, pt1, pt2, pt3, &pt_useless)!=0) continue;
	if(seg_cross_chck(receiv, pt1, pt0, pt2, &pt_useless)==1) continue; 

	if(!eq_pts (pt0, pt1))
	    if(seg_cross_chck(receiv, pt2, pt0, pt1, &pt_useless)==1) continue;
	if(seg_cross_chck(receiv, pt2, pt3, pt1, &pt_useless)!=0) continue; 

        if(!eq_pts (pt0, pt1))
	    if(seg_cross_chck(receiv, pt3, pt0, pt1, &pt_useless)==1) continue;
	if(seg_cross_chck(receiv, pt3, pt2, pt0, &pt_useless)==1) continue; 

	rec[sp++] = rec[ii];
	out->flag[rec[ii]] ++;
    }
    *num_rec = sp;
    return;
}

/*----------------------------------------------------------------------------*/

static int seg_cross_chck (struct point pt0, struct point pt1, struct point pp0, struct point pp1, 
			   struct point *pt)
/* Checks if segment pt0-pt1 crosses (intersects) segment pp0-pp1.
* The subroutine returns a one in case they intersect and zero
* in case they don't. */
{
    float mt=0.0f, bt=0.0f, mp=0.0f, bp=0.0f;
    float X, Z;
    int inft, infp;

    inft = infp = 0;

    /* Obtain line eqn. that defines both segments */
    if(pt0.x == pt1.x) inft = 1;
    else {
	mt = slope (pt0, pt1);
	bt = ordinate (pt0, mt);
    }

    if(pp0.x == pp1.x) infp = 1;

    else {
	mp = slope (pp0, pp1);
	bp = ordinate (pp0, mp);
    }

    /* Obtain intersection point of the two lines */

    if(inft && infp) return 0;
    if(inft) {
	X = pt0.x;
	Z = mp * X + bp;
    } else if(infp) {
	X = pp0.x;
	Z = mt * X + bt;
    } else { 
	if (mt==mp) return 0;
        X = (bt - bp) / (mp - mt);
        Z = mt * X + bt;
    }
    *pt = makepoint (X, Z);

    /* Check if it belongs to both segments */

    if(!bel(X, pt0.x, pt1.x)) return 0;
    if(!bel(Z, pt0.z, pt1.z)) return 0;
    if(!bel(X, pp0.x, pp1.x)) return 0;
    if(!bel(Z, pp0.z, pp1.z)) return 0;

    /* If pt0 (receiver) is equal to the intersection point: */
    if(eq_pts(pt0, *pt)) { 
	if((pt0.x<=gox) || (pt0.x>=gox+(gnx-1)*gdx) ||                                     (pt0.z<=goz) || (pt0.z>=goz+(gnz-1)*gdz)) return 0;
	return 2;
    }

    /* At this point we know that both segments intersect */
    return 1;
}

/*----------------------------------------------------------------------------*/

static void interp_rec(int *rec, struct grid *out, int num_rec, 
		       struct heptagon c0, struct heptagon c1, struct heptagon c2, struct heptagon c3, 
		       float tim, float *vel, int first)
{
    int ii;
    struct point pt0, pt1, pt2, pt3;
    struct point ptr, wfront_ptr;
    float d1, d2, s, vmed, ampl;
    float d0r, d2r, d02;
    float d1r, d3r, d13;

    pt0 = c1.x0; pt1 = c2.x0;
    pt2 = c1.x1; pt3 = c2.x1;
    for(ii=0;ii<num_rec;ii++) {
	ptr = makepoint(ROUND(rec[ii]/gnz)*gdx+gox, (rec[ii]%gnz)*gdz+goz);

	d0r = dist(pt0, ptr);
	d2r = dist(pt2, ptr);
	d02 = dist(pt0, pt2);

	d1 = (d0r*d0r + d02*d02 - d2r*d2r) / (2. * d02);
	d1 = sqrt(fabs((double) d0r*d0r - d1*d1));

	d1r = dist(pt1, ptr);
	d3r = dist(pt3, ptr);
	d13 = dist(pt1, pt3);

	d2 = (d1r*d1r + d13*d13 - d3r*d3r) / (2. * d13);
	d2 = sqrt(fabs((double) d1r*d1r - d2*d2));

	s = d1 / (d1 + d2);

	wfront_ptr = ptinterp (c0.x0, pt0, pt1, c3.x0, s);
	vmed =  velo (ptr.x, ptr.z, vel);
	vmed += velo (wfront_ptr.x ,wfront_ptr.z , vel);

	ampl =realinterp(c0.x0,pt0,pt1,c3.x0,c0.ampl,c1.ampl,c2.ampl,c3.ampl,s);
	ampl *= sqrt(dist(pt0, pt1)/ (d1+d2));

	// Added direction cosines - note there might be better ways to do this
	float angle = realinterp(c0.x0,pt0,pt1,c3.x0,c0.angle,c1.angle,c2.angle,c3.angle,s);
	float dirx = sin(angle);
	float dirz = cos(angle);

	// Direction cosines at source
	float srcAng = realinterp(c0.x0,pt0,pt1,c3.x0,c0.srcAng,c1.srcAng,c2.srcAng,c3.srcAng,s);
	float srcx = sin(srcAng);
	float srcz = cos(srcAng);


	// The amplitude term is just computin ray path length
	//
	// Compute Proper geometrical spreading - dOmega/dA where Omega = angle at source, dA - arc length
	//

	// we have angles in range 0-2pi
	float minth = SF_MIN(c1.srcAng,c2.srcAng);
	float maxth = SF_MAX(c1.srcAng,c2.srcAng);
	float dt1 =  maxth-minth; 							// standard angle difference
	float dt2 =  2.*SF_PI-dt1;							// angle difference assuming wrap around
	float dTheta = SF_MIN( dt1,  dt2 );					// Shortest angle difference between the two rays
	float dArc = d1+d2;									// arc length at receiver

	float offset = fabs(wfront_ptr.x-sourcePosition.x); // Offset from the source position

	// Compute wavefront solid angle at source
	float dOmega;
	float dArea;
	if((SF_PI >= minth && SF_PI <= maxth)  || dTheta==dt2){
		// If vertical  - treat as cap
		dOmega=2.*SF_PI*(1.-cos(dTheta/2.));
		dArea=SF_PI*(dArc/2.)*(dArc/2.);
	}else{
		// treat as rectangle going around the sphere.
		float c1 = minth;
		if(minth>SF_PI) c1 = 2.*SF_PI - minth;
		float c2 = maxth;
		if(maxth>SF_PI) c2 = 2.*SF_PI - maxth;
		float colatN = SF_MIN(c1,c2);
		float colatS =  SF_MAX(c1,c2);
		dOmega=(cos(colatN)-cos(colatS))*2*SF_PI;
		dArea = dArc*2*SF_PI*offset;
	}

	// We now save (Note - in theory this is = R = ampl in a homogenous medium)
	float invGeoSpread=sqrt(dArea/dOmega);
	//float invGeoSpread=dArea/dOmega;

	if(first) {
	    if(out->flag[rec[ii]] > 1) {
	    	float told=out->time[rec[ii]];
	    	out->time[rec[ii]] = SF_MIN(out->time[rec[ii]],2.* dist(ptr,wfront_ptr) / vmed + tim);
	    	out->ampl[rec[ii]] = SF_MIN(out->ampl[rec[ii]], ampl);

	    	if(out->time[rec[ii]] != told){
	    		//out->ampl[rec[ii]] = ampl;
	    		out->dirx[rec[ii]] = dirx;
		    	out->dirz[rec[ii]] = dirz;
	    		out->srcx[rec[ii]] = srcx;
		    	out->srcz[rec[ii]] = srcz;
		    	out->invgeo[rec[ii]] = invGeoSpread;
	    	}
		} else {
	    	out->time[rec[ii]] = 2.* dist(ptr,wfront_ptr) / vmed + tim;
	    	out->ampl[rec[ii]] = ampl;
		    out->dirx[rec[ii]] = dirx;
		    out->dirz[rec[ii]] = dirz;
		    out->srcx[rec[ii]] = srcx;
		    out->srcz[rec[ii]] = srcz;
	    	out->invgeo[rec[ii]] = invGeoSpread;
		}
	} else {
	    out->time[rec[ii]] = 2. * dist(ptr, wfront_ptr) / vmed + tim;
	    out->ampl[rec[ii]] = ampl;
	    out->dirx[rec[ii]] = dirx;
	    out->dirz[rec[ii]] = dirz;
	    out->srcx[rec[ii]] = srcx;
	    out->srcz[rec[ii]] = srcz;
    	out->invgeo[rec[ii]] = invGeoSpread;
	}
    }
    return;
}

/*----------------------------------------------------------------------------*/
