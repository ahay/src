/* Subroutines for VR programs */
/*
  Copyright (C) 1991 The Board of Trustees of Stanford University
  
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

#include "vrsubs.h"
#include "rot.h"

/*
 * How much to perturb when taking numerical derivatives.
 * Shouldn't be too small compared to machine precision,
 * or you're in trouble! On the other hand, the smallest allowed
 * tile size in the slowness domain is limited by how big it is.
 * ("Hit DIS limit" error message.)
 */
#define DIS  (.00000001)/* Derivative interval, in degrees */

#ifndef _vrsubs_h

struct wave
{
/* Normalized particle motion direction ({X, Y, Z} components) */
    double          part[3];
/*
 *  Phase Slowness 
 */
    double          pslo;
};
/*^*/

struct imp
{
/* Normalized X, Y, and Z coordinates in Group domain */
    double          co[3];
/* Group Slowness */
    double          gslo;
};
/*^*/

struct trip
{
/* Normalized X, Y, and Z coordinates in Slowness domain */
    double          co[3];
/* Stuff specific to each wavetype */
    struct wave     wv[3];
};
/*^*/

struct corner
{
/*
 * A triplet of points in slowness domain; 0 is the actual corner,
 * 1 and 2 are only used for taking numerical derivatives
 */
    struct trip     tr[3];
/*
 * Normalized coordinate and slowness of point corresponding to
 * given corner in GROUP domain for each of 3 wavetypes.
 */
    struct imp      im[3];
};
/*^*/

/* Color */
struct color
{
    float           red;
    float           green;
    float           blue;
};
/*^*/

struct patch
{
/*
 * 1-4 are four corners, 0 is center
 * Due to programmer sloppiness, not all parts of cn[0] are actually used!
 */
    struct corner   cn[5];
/* Whether given edge is part of axis or not */
    int             drawmask[4];
/*
 * How many times tile has been subdivided in each of 2 directions
 * (latitude and longitude). Tiles containing singularities will always
 * "need to be subdivided". This allows the program to keep track of how
 * many times it has subdivided already so it can eventually give up.
 */
    int             level[2];
/* Color of particle motion direction for each of 3 wavetypes */
    struct color    color[3];
/* Whether tile is a candidate to be deleted */
    bool             clipit;
};
/*^*/


#endif

/*
 * Matrix of Particle motion directions (eigenvectors)
 * A is which wavetype, B is which component {x,y,z}
 */
#define	PART(A,B) part[3*(A)+(B)]
static    double  part[3 * 3];

void phaseslo (const double *cc, double xx, double yy, double zz, double *pslo)
/*< Given a normalized direction (xx,yy,zz) in Phase Slowness space,
 * return the phase slownesses of the 3 wavetypes. >*/
{
    double pvelret[3];
    int             ii;
    int             error;

/*
 * The actual Matrix you take the eigenvalues and eigenvectors of
 * to find out the 3 wavetypes' velocities and particle motion directions.
 */
#define	CHRIS(A,B)  chris[3*(A)+(B)]
    double          chris[3 * 3];


    error = 0;
    impulse_ (cc, &xx, &yy, &zz, pvelret, chris, part, &error);

    if (error)
    {
	fprintf (stderr, "x=%f y=%f z=%f\n", xx, yy, zz);
	fprintf (stderr, "\n");
	fprintf (stderr, "%.24f %.24f %.24f\n", CHRIS (0, 0), CHRIS (1, 0), CHRIS (2, 0));
	fprintf (stderr, "%.24f %.24f %.24f\n", CHRIS (0, 1), CHRIS (1, 1), CHRIS (2, 1));
	fprintf (stderr, "%.24f %.24f %.24f\n", CHRIS (0, 2), CHRIS (1, 2), CHRIS (2, 2));
	fprintf (stderr, "\n");
	fprintf (stderr, "%f %f %f\n", pvelret[0], pvelret[1], pvelret[2]);
	fprintf (stderr, "\n\n");

	fprintf (stderr, "SVD CRASH %d!\n", error);
	exit (1);
    }

    for (ii = 0; ii < 3; ii++)
    {
	pslo[ii] = 1. / pvelret[ii];
    }
}

void xyz (double ar, double up, double *xx, double *yy, double *zz)
/*< Go from (Lat, Long) to (x,y,z) >*/
{
    *xx = cos (up) * cos (ar);
    *yy = cos (up) * sin (ar);
    *zz = sin (up);
}


void swapit (const double *cc, struct trip *trip, struct imp *imp)
/*< This is the main place where the mathematics occurs.
 * This routine takes a point in phase slowness space, and
 * calculates the corresponding point in group slowness space. >*/
{
    double          x1, y1, z1, pslo1[3];
    double          x2, y2, z2, pslo2[3];

    double          xn, yn, zn, gslo;
    double          x0, y0, z0;

    double          tt0, tt1, norm;
    int             ii, jj;

/*
 * Get the coordinates of the point
 * in question in phase slowness space.
 */
    x0 = trip[0].co[0];
    y0 = trip[0].co[1];
    z0 = trip[0].co[2];

/*
 * Perturb away along a line of longitude to get another point
 */
    tt1 = tt0 = sqrt (x0 * x0 + y0 * y0);
    z1 = z0;
    rot (DIS, &tt1, &z1);
    if (tt0 == 0.)
    {
/*
 * Be careful if original point had Lat = +- 90,
 * just arbitrarily pick to go along X axis in that case.
 */
	x1 = tt1;
	y1 = 0.;
    }
    else
    {
	x1 = x0 * tt1 / tt0;
	y1 = y0 * tt1 / tt0;
    }
/*
 * Calculate Christoffel stuff for this perturbed point.
 */
    phaseslo (cc, x1, y1, z1, pslo1);

/*
 * Save all the info in the 1st slot in the "trip" structure.
 */
/* Coordinate */
    trip[1].co[0] = x1;
    trip[1].co[1] = y1;
    trip[1].co[2] = z1;
/* Loop over 3 wavetypes */
    for (ii = 0; ii < 3; ii++)
    {
/* Phase Slowness */
	trip[1].wv[ii].pslo = pslo1[ii];
	for (jj = 0; jj < 3; jj++)
	{
/* Particle Motion direction */
	    trip[1].wv[ii].part[jj] = PART (ii, jj);
	}
    }

/*
 * We are going to use these auxilliary points to do numerical
 * derivatives. We need to make sure we are comparing points on
 * the same surface when we do so. "orderit" swaps around the
 * ordering of the 3 wavetypes so that the particle motion directions
 * match the information in slot 0 of the "trip" structure.
 */
    orderit (trip[0].wv, trip[1].wv);

/*
 * Coming up with the next auxilliary point is harder.
 * We want it to be perturbed from the original point in a
 * perpendicular direction to the perturbation for the other
 * auxilliary point. We do it with a cross product.
 * Take (auxilliary point - original point) and cross that
 * with (original point - origin) to get perpendicular perturbation,
 * then add that to the original point to get the new point.
 */
    x2 = y0 * (z1 - z0) - z0 * (y1 - y0) + x0;
    y2 = -x0 * (z1 - z0) + z0 * (x1 - x0) + y0;
    z2 = x0 * (y1 - y0) - y0 * (x1 - x0) + z0;
/* Then normalize. */
    norm = sqrt (x2 * x2 + y2 * y2 + z2 * z2);
    x2 /= norm;
    y2 /= norm;
    z2 /= norm;
/* Calculate Christoffel stuff for new auxilliary point */
    phaseslo (cc, x2, y2, z2, pslo2);

/* Save the new information in slot 2 of "trip". */
    trip[2].co[0] = x2;
    trip[2].co[1] = y2;
    trip[2].co[2] = z2;
    for (ii = 0; ii < 3; ii++)
    {
	trip[2].wv[ii].pslo = pslo2[ii];
	for (jj = 0; jj < 3; jj++)
	{
	    trip[2].wv[ii].part[jj] = PART (ii, jj);
	}
    }

/* Again order wavetypes to match original point's ordering. */
    orderit (trip[0].wv, trip[2].wv);

/*
 * Now for each of 3 wavetypes calculate the numerical derivatives
 * and see where the point in question maps to in the group domain.
 */
    for (ii = 0; ii < 3; ii++)
    {
/* Un-normalize everything */
	x0 = trip[0].co[0] * trip[0].wv[ii].pslo;
	y0 = trip[0].co[1] * trip[0].wv[ii].pslo;
	z0 = trip[0].co[2] * trip[0].wv[ii].pslo;

	x1 = trip[1].co[0] * trip[1].wv[ii].pslo;
	y1 = trip[1].co[1] * trip[1].wv[ii].pslo;
	z1 = trip[1].co[2] * trip[1].wv[ii].pslo;

	x2 = trip[2].co[0] * trip[2].wv[ii].pslo;
	y2 = trip[2].co[1] * trip[2].wv[ii].pslo;
	z2 = trip[2].co[2] * trip[2].wv[ii].pslo;

/*
 * Find group direction.
 * Group direction is perpendicular to Phase Slowness surface.
 * We have constructed a little right triangle on this surface;
 * to get perpendicular we take cross product of the two edges
 * of the triangle meeting at the right angle. (pt1 - pt0) X (pt2 - pt0)
 */
	xn = (y1 - y0) * (z2 - z0) - (z1 - z0) * (y2 - y0);
	yn = -(x1 - x0) * (z2 - z0) + (z1 - z0) * (x2 - x0);
	zn = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);
/*
 * Normalize group direction.
 */
	norm = sqrt (xn * xn + yn * yn + zn * zn);
	xn /= norm;
	yn /= norm;
	zn /= norm;
/*
 * Now find Group Slowness.
 * This is the component of the Phase Slowness vector in the
 * group direction. Just take the dot product of the Phase Slowness
 * vector (pt0) with the normalized group direction vector (ptn).
 */
	gslo = xn * x0 + yn * y0 + zn * z0;

/*
 * Save the results in the "imp" structure.
 */
	imp[ii].co[0] = xn;
	imp[ii].co[1] = yn;
	imp[ii].co[2] = zn;

	imp[ii].gslo = gslo;

/* End loop over 3 wavetypes */
    }
}

void dotile (bool order, bool which, const double *cc, 
	     struct corner *cn, struct color color[3])
/*< Given the coordinates of the 4 corners of a Wyoming-shaped tile
 * in Slowness space, fill in everything else. >*/
{
    float f_blue, f_green, f_red;
    int             ii, jj, kk;
    double          x0, y0, z0, norm;
    double          pslo[3];
/*
 * Matrix giving pure "P", "SV", and "SH" particle directions at a given
 * point.
 */
#define	ORTHO(A,B)  ortho[3*(A)+(B)]
    double          ortho[3 * 3];

/*
 * Find center of tile (average of corners 1 through 4),
 * put that (normalized) in corner 0.
 */
    x0 = 0.;
    y0 = 0.;
    z0 = 0.;
    for (jj = 1; jj < 5; jj++)
    {
	x0 += cn[jj].tr[0].co[0];
	y0 += cn[jj].tr[0].co[1];
	z0 += cn[jj].tr[0].co[2];
    }

    norm = sqrt (x0 * x0 + y0 * y0 + z0 * z0);
    x0 /= norm;
    y0 /= norm;
    z0 /= norm;

    cn[0].tr[0].co[0] = x0;
    cn[0].tr[0].co[1] = y0;
    cn[0].tr[0].co[2] = z0;

/*
 * Now we have the direction to the center of the tile.
 * Save it as the "pure P" direction.
 */
    ORTHO (0, 0) = x0;
    ORTHO (0, 1) = y0;
    ORTHO (0, 2) = z0;

/*
 * "Pure SH" direction.
 * The center of a tile should never have |z| = 1,
 * so the division is safe.
 */
    norm = sqrt (x0 * x0 + y0 * y0);
    ORTHO (2, 0) = -y0 / norm;
    ORTHO (2, 1) = x0 / norm;
    ORTHO (2, 2) = 0.;

/*
 * "Pure SV" direction. Cross product of the P and SH directions.
 */
    ORTHO (1, 0) = ORTHO (0, 1) * ORTHO (2, 2) - ORTHO (0, 2) * ORTHO (2, 1);
    ORTHO (1, 1) = -(ORTHO (0, 0) * ORTHO (2, 2) - ORTHO (0, 2) * ORTHO (2, 0));
    ORTHO (1, 2) = ORTHO (0, 0) * ORTHO (2, 1) - ORTHO (0, 1) * ORTHO (2, 0);

/*
 * Loop over center and 4 corners
 */
    for (jj = 0; jj < 5; jj++)
    {
	x0 = cn[jj].tr[0].co[0];
	y0 = cn[jj].tr[0].co[1];
	z0 = cn[jj].tr[0].co[2];

/*
 * Solve Christoffel equation for each; save all the resulting
 * information. (Phase Slowness, particle motion directions)
 */
	phaseslo (cc, x0, y0, z0, pslo);
/* Loop over 3 wavetypes */
	for (ii = 0; ii < 3; ii++)
	{
	    cn[jj].tr[0].wv[ii].pslo = pslo[ii];

/* Loop over {X,Y,Z} */
	    for (kk = 0; kk < 3; kk++)
	    {
		cn[jj].tr[0].wv[ii].part[kk] = PART (ii, kk);
	    }
	}

/* Center is special case, skip it */
	if (jj > 0)
	{
/*
 * If "order=y" option, try to reorder order of wavetypes at
 * each corner to match particle motions with center. This is
 * a sensible thing to do for TI, so you can ignore possible
 * "intersection" singularities. Of course for general anisotropy
 * there is no such thing as an "intersection" singularity. Bad
 * idea then.
 */
	    if (order)
		orderit (cn[0].tr[0].wv, cn[jj].tr[0].wv);

/*
 * Calculate the image of each corner in the group slowness domain.
 * This can be computationally expensive, so only do it if asked to.
 */
	    if (which)
		swapit (cc, cn[jj].tr, cn[jj].im);
	}
    }

/*
 * Calculate the unshaded color for the tile, based on the particle
 * motion direction at the center.
 */
/* Loop over 3 wavetypes */
    for (jj = 0; jj < 3; jj++)
    {
/* P */
	f_blue = fabs (dot (&(ORTHO (0, 0)),
			    cn[0].tr[0].wv[jj].part));
/* SH */
	f_green = fabs (dot (&(ORTHO (1, 0)),
			     cn[0].tr[0].wv[jj].part));
/* SV */
	f_red = fabs (dot (&(ORTHO (2, 0)),
			   cn[0].tr[0].wv[jj].part));

	color[jj].red = f_red;
	color[jj].green = f_green;
	color[jj].blue = f_blue;
    }
}

void subdivide1 (bool order, bool which, const double *cc, 
		 struct patch *pt1, struct patch *pt2)
/*< Split one patch into 2 over longitude. "pt1" is input patch that will be split in
 * two, with new half going into "pt2". >*/
{
    int             ii;
    double          cotemp[3], norm;

/*
 * Copy info across from old patch to new
 */
    (*pt2).clipit = (*pt1).clipit;

    (*pt1).level[0]++;
    (*pt2).level[0] = (*pt1).level[0];
    (*pt2).level[1] = (*pt1).level[1];

    for (ii = 0; ii < 3; ii++)
    {
	(*pt2).cn[1].tr[0].co[ii] =
	    (*pt1).cn[1].tr[0].co[ii];
	(*pt2).cn[4].tr[0].co[ii] =
	    (*pt1).cn[4].tr[0].co[ii];
    }

    (*pt2).drawmask[0] = (*pt1).drawmask[0];
    (*pt2).drawmask[1] = 0;
    (*pt2).drawmask[2] = (*pt1).drawmask[2];
    (*pt2).drawmask[3] = (*pt1).drawmask[3];
    (*pt1).drawmask[3] = 0;

/*
 * Calculate 2 new corners, put them into new patches
 */
    norm = 0.;
    for (ii = 0; ii < 3; ii++)
    {
	cotemp[ii] =
	    ((*pt1).cn[1].tr[0].co[ii]
	     +
	     (*pt1).cn[2].tr[0].co[ii]) / 2.;
	norm += cotemp[ii] * cotemp[ii];
    }
    norm = sqrt (norm);
    for (ii = 0; ii < 3; ii++)
    {
	(*pt1).cn[1].tr[0].co[ii] =
	    (*pt2).cn[2].tr[0].co[ii] =
	    cotemp[ii] / norm;
    }


    norm = 0.;
    for (ii = 0; ii < 3; ii++)
    {
	cotemp[ii] =
	    ((*pt1).cn[4].tr[0].co[ii]
	     +
	     (*pt1).cn[3].tr[0].co[ii]) / 2.;
	norm += cotemp[ii] * cotemp[ii];
    }
    norm = sqrt (norm);
    for (ii = 0; ii < 3; ii++)
    {
	(*pt1).cn[4].tr[0].co[ii] =
	    (*pt2).cn[3].tr[0].co[ii] =
	    cotemp[ii] / norm;
    }

/*
 * Recalculate everything 2 new modified patches.
 * This repeats some computations unnecissarily at the 4 corners
 * that stayed the same, but tough. It's not worth the extra
 * code complexity to keep track of what's been calculated
 * and what hasn't. (Remember you also have the additional
 * complication of possibly having to reorder _everything_
 * to match the new center's wavetype order, too. Ugh.)
 */
    dotile (order, which, cc, (*pt1).cn, (*pt1).color);
    dotile (order, which, cc, (*pt2).cn, (*pt2).color);
}

void subdivide2 (bool order, bool which, const double *cc,
		 struct patch *pt1, struct patch *pt2)
/*< Just like subdivide1, but splitting over latitude. >*/
{
    int             ii;
    double          cotemp[3], norm;

    (*pt2).clipit = (*pt1).clipit;

    (*pt1).level[1]++;
    (*pt2).level[1] = (*pt1).level[1];
    (*pt2).level[0] = (*pt1).level[0];

    for (ii = 0; ii < 3; ii++)
    {
	(*pt2).cn[4].tr[0].co[ii] =
	    (*pt1).cn[4].tr[0].co[ii];
	(*pt2).cn[3].tr[0].co[ii] =
	    (*pt1).cn[3].tr[0].co[ii];
    }

    (*pt2).drawmask[0] = 0;
    (*pt2).drawmask[1] = (*pt1).drawmask[1];
    (*pt2).drawmask[2] = (*pt1).drawmask[2];
    (*pt2).drawmask[3] = (*pt1).drawmask[3];
    (*pt1).drawmask[2] = 0;

    norm = 0.;
    for (ii = 0; ii < 3; ii++)
    {
	cotemp[ii] =
	    ((*pt1).cn[1].tr[0].co[ii]
	     +
	     (*pt1).cn[4].tr[0].co[ii]) / 2.;
	norm += cotemp[ii] * cotemp[ii];
    }
    norm = sqrt (norm);
    for (ii = 0; ii < 3; ii++)
    {
	(*pt1).cn[4].tr[0].co[ii] =
	    (*pt2).cn[1].tr[0].co[ii] =
	    cotemp[ii] / norm;
    }


    norm = 0.;
    for (ii = 0; ii < 3; ii++)
    {
	cotemp[ii] =
	    ((*pt1).cn[2].tr[0].co[ii]
	     +
	     (*pt1).cn[3].tr[0].co[ii]) / 2.;
	norm += cotemp[ii] * cotemp[ii];
    }
    norm = sqrt (norm);
    for (ii = 0; ii < 3; ii++)
    {
	(*pt1).cn[3].tr[0].co[ii] =
	    (*pt2).cn[2].tr[0].co[ii] =
	    cotemp[ii] / norm;
    }

    dotile (order, which, cc, (*pt1).cn, (*pt1).color);
    dotile (order, which, cc, (*pt2).cn, (*pt2).color);
}


void orderit (struct wave *waveref, struct wave *waveshf)
/*< Reorder order of wavetypes in "waveshf" to make closest
 * match to wavetypes in "waveref". >*/
{
    struct wave     temp[3];
    int             ii, jj, ll;
    double          sum[3];
    int             which;
    double          val;

/*
 * For each of three reference wavetypes, find the one with the
 * closest particle motion direction.
 */
/* Loop over reference particle motion direction */
    for (ll = 0; ll < 3; ll++)
    {
/* Loop to find best match to it */
	for (ii = 0; ii < 3; ii++)
	{
	    sum[ii] = 0.;
	    for (jj = 0; jj < 3; jj++)
	    {
		sum[ii] += waveref[ll].part[jj] * waveshf[ii].part[jj];
	    }
	    sum[ii] = fabs (sum[ii]);
	}

	which = 0;
	val = sum[0];
	for (ii = 1; ii < 3; ii++)
	{
	    if (sum[ii] > val)
	    {
		which = ii;
		val = sum[ii];
	    }
	}
	temp[ll] = waveshf[which];
    }

/* Copy best matches across */
    for (ll = 0; ll < 3; ll++)
    {
	waveshf[ll] = temp[ll];
    }
}
