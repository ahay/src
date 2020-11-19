/* Plot impulse responses in 3 dimensions */
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
#include <stdlib.h>

#include <rsf.h>
#include <rsfplot.h>

#include "vrsubs.h"
#include "rot.h"

#define TOPI (2.*SF_PI)
#define PITO (.5*SF_PI)

/*
 * How big to initially allocate space for tiles
 */
#define NPOINTS 1000	/* "work" (patch) tiles, containing everything */
#define NPATCH (NPOINTS*2)	/* stripped-down "graphical output"
				 * (polypatch) tiles */

/*
 * How much to perturb when taking numerical derivatives.
 * Shouldn't be too small compared to machine precision,
 * or you're in trouble! On the other hand, the smallest allowed
 * tile size in the slowness domain is limited by how big it is.
 * ("Hit DIS limit" error message.)
 */
#define DIS  (.00000001)/* Derivative interval, in degrees */

/* Christoffel Matrix */
#define	CC(A,B)	cc[6*(A)+(B)]
static double   cc[6 * 6];

/*
 * This is the structure that contains all the information needed to
 * plot a 4-sided (Wyoming-shaped) tile ("patch") of the slowness surface
 * and the corresponding piece of the impulse surface.
 */
struct polypatch
{
/*
 * x, y, and z coordinates of the Phase Slowness surface.
 * 1-4 are the 4 corners; 0 is the center (average of the other 4).
 * The 4 corners are normalized (you need to multiply by "pslo"),
 * the center is not.
 */
    double          xx[5], yy[5], zz[5];
/*
 * pslo[1-5] gives the Phase Slowness for the 4 corners.
 * pslo[0] is YES if the Phase Slowness information is present,
 * NO otherwise.
 */
    double          pslo[5];
/*
 * Ditto for the Group Velocity surface. Note the 4 corners on the
 * two surfaces correspond, but the two centers need not.
 */
    double          xxn[5], yyn[5], zzn[5];
    double          gvel[5];
/*
 * The particle motion direction at the center of the slowness surface
 * tile.
 */
    double          part[3];
/*
 * If YES, then this tile has particle motion directions at at least
 * one of the corners that deviates more than "coslimit" from the
 * particle motion at the center. This almost certainly means this
 * tile contains a singularity.
 */
    int             singular;
/*
 * If 2, then the given edge of the tile is one of the axis lines
 * that usually is drawn white and fat. If 1, then an original
 * tile edge. 0 if created by later subdivision.
 */
    int             drawmask[4];
/*
 * Which surface is this tile on, 0=fastest 2=slowest
 */
    int             kind;
/*
 * Color of this tile before shading
 * (red=SH, green=SV, blue=P)
 */
    float           red;
    float           green;
    float           blue;
/*
 * Latitude and Longitude of center of this tile in Slowness domain.
 * Latitude +90 = Z axis; Longitude 0 = X axis; Longitude 90 = Y axis
 */
    float           lati;
    float           longi;
};

static int pt_len;
static struct patch *pt;

static int points_len;
static struct polypatch *points;

static void morepatches(void);
static void morepolys(void);

int main(int argc, char* argv[])
{
    double          cosdis;	/* Lower limit is set to 10 times DIS */

    float           start_angle, end_angle, bottom_angle, top_angle;
    bool            which;
    bool            order;

    float           zmin, zmax;
    float           ymin, ymax;
    float           xmin, xmax;

    bool            or;
    int             skip;

    float           zzmin, zzmax;
    float           yymin, yymax;
    float           xxmin, xxmax;

    int             what[3];
    int             iinc;
    int             iinc2;

    int             max_level;
    float           coslimit;

    char    *singstring;
    FILE    *singout;

    int             i, count, ncount;
    double          ar, up, up2, incup, cosincup2, cosincup2x2, incar, tinc;
    double          ara[5], upa[5];
    double          x, y, z;
    double          xav, yav, zav;
    double          xnav, ynav, znav;
    double          x0, y0, z0, pslo0;
    double          xn, yn, zn, gslo;
    char            string[80];
    int             npoints, ii, jj, discarded;
    float           temp, ccnorm;
    double          d1, d2;
    bool            dislim1, dislim2;
    bool            subit;
    double          sum;
    bool            subdivided;
    double          cosangle, cosangle1, cosangle2;

    int patchlen;
    unsigned char *raw;
    sf_file out;

    sf_init(argc,argv);
    out = sf_output("out");

    pt_len = NPOINTS;
    pt = (struct patch *) sf_alloc (pt_len, sizeof (struct patch));

    points_len = NPATCH;
    points = (struct polypatch *) sf_alloc (points_len,sizeof (struct polypatch));

    if (NULL == (singstring = sf_getstring ("sing"))) singstring="Vr3d.log";
    /* Log file */

    singout = fopen (singstring, "w");
    if (singout == NULL) sf_error( "Can't open log file.");
 
    if (!sf_getfloat ("start",&start_angle)) start_angle = 0.;
    /* longitude start */

    if (!sf_getfloat ("end", &end_angle)) end_angle = 360.;
    /* longitude end */

    if (!sf_getfloat ("bottom", &bottom_angle))  bottom_angle = -90.;
    /* latitude start */
    
    if (!sf_getfloat ("top", &top_angle)) top_angle = 90.;
    /* latitude end */

    if (!sf_getbool ("order", &order)) order=false;
    /* try to swap around the surfaces to make them continuous */

    if (!sf_getint ("maxlevel", &max_level)) max_level=5;
    /* maximum number of re-subdivisions */

    if (!sf_getfloat ("coslimit", &coslimit)) coslimit=25.;
    /* aximum deviation in particle motion angle */

    if (!sf_getfloat ("xmin", &xmin)) xmin=-100.;
    if (!sf_getfloat ("ymin", &ymin)) ymin=-100.;
    if (!sf_getfloat ("zmin", &zmin)) zmin=-100.;
    if (!sf_getfloat ("xmax", &xmax)) xmax=100.;
    if (!sf_getfloat ("ymax", &ymax)) ymax=100.;
    if (!sf_getfloat ("zmax", &zmax)) zmax=100.;


    if (!sf_getbool ("or", &or)) or=false;
    /* modifier: if or=y ORs instead of ANDS the clips. */
 
    if (!sf_getint ("skip", &skip)) skip=-1;
    /*  modifier: skip=-1 don't clip this surface (-1 for none skipped)
     			0 = fastest
     			1 = intermediate
     			2 = slowest
     			3 = red (SH)
     			4 = green (SV)
     			5 = blue (P) */

    if (!sf_getfloat ("xxmin", &xxmin)) xxmin=-100.;
    if (!sf_getfloat ("yymin", &yymin)) yymin=-100.;
    if (!sf_getfloat ("zzmin", &zzmin)) zzmin=-100.;
    if (!sf_getfloat ("xxmax", &xxmax)) xxmax=100.;
    if (!sf_getfloat ("yymax", &yymax)) yymax=100.;
    if (!sf_getfloat ("zzmax", &zzmax)) zzmax=100.;

    if (!sf_getint ("inc", &iinc)) iinc=4;
    /* density of gridding (How many tiles to cover 90 degree of longitude in initial tiling.) */
    if (!sf_getint ("inc2", &iinc2)) iinc2 = iinc;
    /* tiles bigger than 90 deg / iinc2 in any dimension will be subdivided to fit */
    if (!sf_getints ("what", what, 3)) { what[0]=1; what[1]=2; what[2]=2; }
    /* which surfaces to do, fastest to slowest */
    if (!sf_getbool ("which", &which)) which=true;
    /* if y, plot impulse response; if n, plot slowness surface */

    sf_putfloat(out,"xmin", xmin);
    sf_putfloat(out,"ymin", ymin);
    sf_putfloat(out,"zmin", zmin);
    sf_putfloat(out,"xmax", xmax);
    sf_putfloat(out,"ymax", ymax);
    sf_putfloat(out,"zmax", zmax);
    sf_putfloat(out,"xxmin", xxmin);
    sf_putfloat(out,"yymin", yymin);
    sf_putfloat(out,"zzmin", zzmin);
    sf_putfloat(out,"xxmax", xxmax);
    sf_putfloat(out,"yymax", yymax);
    sf_putfloat(out,"zzmax", zzmax);
    sf_putfloat(out,"start", start_angle);
    sf_putfloat(out,"end", end_angle);
    sf_putfloat(out,"bottom", bottom_angle);
    sf_putfloat(out,"top", top_angle);
    sf_putint(out,"or", (or != false));
    sf_putint(out,"skip", skip);

    sf_putint(out,"inc", iinc);
    sf_putint(out,"inc2", iinc2);
    sf_putints(out,"what", what, 3);
    sf_putfloat(out,"coslimit",coslimit);
    coslimit = cosf ( coslimit * TOPI / 360.);
    cosdis = cosf (DIS * 10. * TOPI / 360.);
    sf_putint(out,"which", (which != false));
    sf_putint(out,"order", order);
    sf_putint(out,"maxlevel", max_level);

    if (!sf_getfloat ("norm", &ccnorm)) ccnorm=1.;
    /* amount to divide everything by */
    sf_putfloat(out,"norm", ccnorm);

    for (ii = 0; ii < 6; ii++)
	for (jj = 0; jj <= ii; jj++)
	{
	    snprintf (string, 80, "c%d%d c%d%d", ii + 1, jj + 1, jj + 1, ii + 1);

	    if (!sf_getfloat (string, &temp)) temp=0.0;
	    sf_putfloat (out, string, temp);
	    
	    CC (ii, jj) = CC (jj, ii) = temp / ccnorm;
	}

    fprintf (singout, "#\n#\n#\tStiffness (C) matrix:\n#\n");
    for (ii = 0; ii < 6; ii++)
    {
	fprintf (singout, "#\t");
	for (jj = 0; jj < 6; jj++)
	{
	    fprintf (singout, "%g\t", ccnorm * CC (ii, jj));
	}
	fprintf (singout, "\n");
    }
    fprintf (singout, "#\n#\n");

/*
 * Total number of tiles
 */
    count = 0;

    incup = PITO / iinc;
    if (iinc2 <= 0)
	cosincup2 = 0.;
    else
	cosincup2 = cos (PITO / iinc2);

    cosincup2x2 = 2. * cosincup2 * cosincup2 - 1.;

    start_angle *= TOPI / 360.;
    end_angle *= TOPI / 360.;

    bottom_angle *= TOPI / 360.;
    top_angle *= TOPI / 360.;

/* Latitude loop */
    for (up = bottom_angle + .5 * incup; up < top_angle - .49 * incup;
	 up += incup)
    {
	sf_warning ("%f, %d;", (float) (up * 360. / TOPI), count);

	if (up > 0.)
	    up2 = up - .5 * incup;
	else
	    if (up < 0.)
		up2 = up + .5 * incup;
	    else
		up2 = up;

	tinc = (int) (.5 + 4. * iinc * fabs (cos (up2)));
	incar = 8;
	while (incar < tinc)
	{
	    incar *= 2;
	}
	incar = TOPI / incar;

/* Longitude loop */
	for (ar = start_angle + .5 * incar; ar < end_angle + .5 * incar - 1e-6; ar += incar)
	{
	    ara[0] = ar;
	    upa[0] = up;

	    xyz (ara[0], upa[0],
		 &x, &y, &z);

	    pt[count].clipit = false;

	    if (or)
	    {
		if (
		    (x > xmax || x < xmin) &&
		    (y > ymax || y < ymin) &&
		    (z > zmax || z < zmin)
		    )
		    pt[count].clipit = true;
	    }
	    else
	    {
		if (
		    x > xmax || x < xmin ||
		    y > ymax || y < ymin ||
		    z > zmax || z < zmin
		    )
		    pt[count].clipit = true;
	    }

	    if (pt[count].clipit && skip == -1)
		continue;

/*
 * Wyoming-shaped tile shape and order:
 *
 *	4 ----- 1
 *	|	|
 *	|	|
 *	|	|
 *	3 ----- 2
 */
	    ara[1] = ar + .5 * incar;
	    upa[1] = up + .5 * incup;

	    ara[2] = ar + .5 * incar;
	    upa[2] = up - .5 * incup;

	    ara[3] = ar - .5 * incar;
	    upa[3] = up - .5 * incup;

	    ara[4] = ar - .5 * incar;
	    upa[4] = up + .5 * incup;

/* Save coordinates */
	    for (jj = 1; jj < 5; jj++)
	    {
		xyz (ara[jj], upa[jj],
		     &x0, &y0, &z0);

		pt[count].cn[jj].tr[0].co[0] = x0;
		pt[count].cn[jj].tr[0].co[1] = y0;
		pt[count].cn[jj].tr[0].co[2] = z0;
	    }

/*
 * Calculate everything given only the coordinates
 * in the Slowness domain
 */
	    dotile (order, which, cc, pt[count].cn, pt[count].color);

	    pt[count].level[0] = 0;
	    pt[count].level[1] = 0;

	    for (jj = 0; jj < 4; jj++)
	    {
		pt[count].drawmask[jj] = 1;
	    }

/*
 * See where Axes intersect tiles
 */
	    if (fabs (pt[count].cn[1].tr[0].co[0]) < 1e-5 &&
		fabs (pt[count].cn[2].tr[0].co[0]) < 1e-5)
		pt[count].drawmask[0] = 2;
	    if (fabs (pt[count].cn[3].tr[0].co[0]) < 1e-5 &&
		fabs (pt[count].cn[4].tr[0].co[0]) < 1e-5)
		pt[count].drawmask[2] = 2;

	    if (fabs (pt[count].cn[1].tr[0].co[1]) < 1e-5 &&
		fabs (pt[count].cn[2].tr[0].co[1]) < 1e-5)
		pt[count].drawmask[0] = 2;
	    if (fabs (pt[count].cn[3].tr[0].co[1]) < 1e-5 &&
		fabs (pt[count].cn[4].tr[0].co[1]) < 1e-5)
		pt[count].drawmask[2] = 2;

	    if (fabs (pt[count].cn[2].tr[0].co[2]) < 1e-5)
		pt[count].drawmask[1] = 2;
	    if (fabs (pt[count].cn[4].tr[0].co[2]) < 1e-5)
		pt[count].drawmask[3] = 2;

	    count++;
	    if (count >= pt_len)
	    {
		morepatches ();
	    }
	}
    } 
    sf_warning(".");

/*
 * Refinement stage
 */
    ncount = count;
/*
 * We work through the data keeping 2 active sites at the same time.
 * We start out with one pointer at the beginning, and one pointer
 * just off the end. Whenever we find a tile at the lower pointer (count)
 * that needs subdividing, we leave one of the created halves there
 * and put the other created half at the end (ncount). Then we increment
 * ncount and repeat. The lower pointer is only incremented when we
 * find a tile that needs no further subdivision.
 */

    sf_warning ("%d tiles", ncount);

    if (which)
    {
/*
 * If "which=1",
 * look at size of tile in group domain.
 */
	for (count = 0; count < ncount;)
	{
	    subdivided = false;

/*
 * First check if we need to cut along Latitude.
 */
	    if (pt[count].level[0] < max_level)
	    {
/*
 * Find length of two parallel sides in Slowness domain
 */
		d1 = dot (pt[count].cn[1].tr[0].co,
			  pt[count].cn[2].tr[0].co);

		d2 = dot (pt[count].cn[4].tr[0].co,
			  pt[count].cn[3].tr[0].co);

/*
 * Are the sides too short compared to the length used when doing
 * numerical derivatives?
 */
		if (d1 > cosdis)
		    dislim1 = true;
		else
		    dislim1 = false;
		if (d2 > cosdis)
		    dislim2 = true;
		else
		    dislim2 = false;

/*
 * Default action is to leave it alone.
 */
		subit = false;

/* Loop over 3 wavetypes */
		for (ii = 0; ii < 3; ii++)
		{
/*
 * Skip ones that we don't care about.
 */
		    if (what[ii] == 0)
			continue;

/*
 * Find length of sides of tile in Group domain
 */
		    d1 = dot (pt[count].cn[1].im[ii].co,
			      pt[count].cn[2].im[ii].co);

		    d2 = dot (pt[count].cn[4].im[ii].co,
			      pt[count].cn[3].im[ii].co);

/*
 * Find cosine of angle between particle motion directions at
 * adjacent corners.
 */
		    cosangle1 = fabs (dot (pt[count].cn[1].tr[0].wv[ii].part,
					   pt[count].cn[2].tr[0].wv[ii].part));

		    cosangle2 = fabs (dot (pt[count].cn[4].tr[0].wv[ii].part,
					   pt[count].cn[3].tr[0].wv[ii].part));

/*
 * If tile is too large in group domain, or particle motion direction is
 * too different at opposite ends, subdivide.
 */
/* Check right side of tile */
		    if (d1 < cosincup2 || cosangle1 < coslimit)
		    {
/*
 * Refuse to subdivide further if the tile is already too small
 * compared to "DIS", the numerical derivative perturbation distance.
 */
			if (dislim1)
			{
			    sf_warning ("Hit DIS limit at level %d", pt[count].level[0]);
			}
			else
			{
			    subit = true;
			}
			break;
		    }

/* Ditto for left side of tile */
		    if (d2 < cosincup2 || cosangle2 < coslimit)
		    {
			if (dislim2)
			{
			    sf_warning ("Hit DIS limit at level %d", pt[count].level[0]);
			}
			else
			{
			    subit = true;
			}
			break;
		    }
		}

/*
 * Finally, subdivide the thing if it needs it.
 */
		if (subit)
		{
		    if (ncount >= pt_len)
		    {
			morepatches ();
		    }
		    subdivide1 (order, which, cc, &pt[count], &pt[ncount]);
		    ncount++;
		    subdivided = true;
		    if (ncount % 100 == 0)
		    {
			sf_warning ("%d tiles, %f percent", ncount, 100. * (float) count / (float) ncount);
		    }
		}
	    }

/*
 * Now check if we need to cut along Longitude.
 * Everything is pretty much the same as before...
 */
	    if (pt[count].level[1] < max_level)
	    {
		d1 = dot (pt[count].cn[4].tr[0].co,
			  pt[count].cn[1].tr[0].co);

		d2 = dot (pt[count].cn[3].tr[0].co,
			  pt[count].cn[2].tr[0].co);

		if (d1 > cosdis)
		    dislim1 = true;
		else
		    dislim1 = false;
		if (d2 > cosdis)
		    dislim2 = true;
		else
		    dislim2 = false;

/*
 * Default action is to leave it alone.
 */
		subit = false;
/* Loop over 3 wavetypes */
		for (ii = 0; ii < 3; ii++)
		{
/*
 * Skip ones that we don't care about.
 */
		    if (what[ii] == 0)
			continue;

		    d1 = dot (pt[count].cn[4].im[ii].co,
			      pt[count].cn[1].im[ii].co);

		    d2 = dot (pt[count].cn[3].im[ii].co,
			      pt[count].cn[2].im[ii].co);

		    cosangle1 = fabs (dot (pt[count].cn[4].tr[0].wv[ii].part,
					   pt[count].cn[1].tr[0].wv[ii].part));

		    cosangle2 = fabs (dot (pt[count].cn[3].tr[0].wv[ii].part,
					   pt[count].cn[2].tr[0].wv[ii].part));

		    if (d1 < cosincup2 || cosangle1 < coslimit)
		    {
			if (dislim1)
			{
			    sf_warning ("Hit DIS limit at level %d", pt[count].level[1]);
			}
			else
			{
			    subit = true;
			}
			break;
		    }

		    if (d2 < cosincup2 || cosangle2 < coslimit)
		    {
			if (dislim2)
			{
			    sf_warning ("Hit DIS limit at level %d", pt[count].level[1]);
			}
			else
			{
			    subit = true;
			}
			break;
		    }
		}

		if (subit)
		{
		    if (ncount >= pt_len)
		    {
			morepatches ();
		    }
		    subdivide2 (order, which, cc, &pt[count], &pt[ncount]);
		    ncount++;
		    subdivided = true;
		    if (ncount % 100 == 0)
		    {
			sf_warning ("%d tiles, %f percent", ncount, 100. * (float) count / (float) ncount);
		    }
		}
	    }

	    if (!subdivided)
	    {
		count++;
	    }
	}
	sf_warning ("%d tiles after refinement", ncount);
    }
    else
    {
/*
 * Refinement stage if we only need to worry about phase domain.
 * Same as above, except we look at size of tile in phase domain.
 */
	for (count = 0; count < ncount;)
	{
	    subdivided = false;

	    if (pt[count].level[0] < max_level)
	    {
		d1 = dot (pt[count].cn[1].tr[0].co,
			  pt[count].cn[2].tr[0].co);

		d2 = dot (pt[count].cn[4].tr[0].co,
			  pt[count].cn[3].tr[0].co);

		if (d1 > cosdis)
		    dislim1 = true;
		else
		    dislim1 = false;
		if (d2 > cosdis)
		    dislim2 = true;
		else
		    dislim2 = false;

/*
 * Default action is to leave it alone.
 */
		subit = false;

		d1 = dot (pt[count].cn[1].tr[0].co,
			  pt[count].cn[2].tr[0].co);

		d2 = dot (pt[count].cn[4].tr[0].co,
			  pt[count].cn[3].tr[0].co);

/* Loop over 3 wavetypes */
		for (ii = 0; ii < 3; ii++)
		{
/*
 * Skip ones that we don't care about.
 */
		    if (what[ii] == 0)
			continue;

		    cosangle1 = fabs (dot (pt[count].cn[1].tr[0].wv[ii].part,
					   pt[count].cn[2].tr[0].wv[ii].part));

		    cosangle2 = fabs (dot (pt[count].cn[4].tr[0].wv[ii].part,
					   pt[count].cn[3].tr[0].wv[ii].part));

		    if (d1 < cosincup2 || cosangle1 < coslimit)
		    {
			if (dislim1)
			{
			    sf_warning ("Hit DIS limit at level %d", pt[count].level[0]);
			}
			else
			{
			    subit = true;
			}
			break;
		    }

		    if (d2 < cosincup2 || cosangle2 < coslimit)
		    {
			if (dislim2)
			{
			    sf_warning ("Hit DIS limit at level %d", pt[count].level[0]);
			}
			else
			{
			    subit = true;
			}
			break;
		    }
		}

		if (subit)
		{
		    if (ncount >= pt_len)
		    {
			morepatches ();
		    }
		    subdivide1 (order, which, cc, &pt[count], &pt[ncount]);
		    ncount++;
		    subdivided = true;
		    if (ncount % 100 == 0)
		    {
			sf_warning ("%d tiles, %f percent", ncount, 100. * (float) count / (float) ncount);
		    }
		}
	    }

	    if (pt[count].level[1] < max_level)
	    {
		d1 = dot (pt[count].cn[4].tr[0].co,
			  pt[count].cn[1].tr[0].co);

		d2 = dot (pt[count].cn[3].tr[0].co,
			  pt[count].cn[2].tr[0].co);

		if (d1 > cosdis)
		    dislim1 = true;
		else
		    dislim1 = false;
		if (d2 > cosdis)
		    dislim2 = true;
		else
		    dislim2 = false;

/*
 * Default action is to leave it alone.
 */
		subit = false;

		d1 = dot (pt[count].cn[4].tr[0].co,
			  pt[count].cn[1].tr[0].co);

		d2 = dot (pt[count].cn[3].tr[0].co,
			  pt[count].cn[2].tr[0].co);

/* Loop over 3 wavetypes */
		for (ii = 0; ii < 3; ii++)
		{
/*
 * Skip ones that we don't care about.
 */
		    if (what[ii] == 0)
			continue;

		    cosangle1 = fabs (dot (pt[count].cn[4].tr[0].wv[ii].part,
					   pt[count].cn[1].tr[0].wv[ii].part));

		    cosangle2 = fabs (dot (pt[count].cn[3].tr[0].wv[ii].part,
					   pt[count].cn[2].tr[0].wv[ii].part));

		    if (d1 < cosincup2 || cosangle1 < coslimit)
		    {
			if (dislim1)
			{
			    sf_warning ("Hit DIS limit at level %d", pt[count].level[1]);
			}
			else
			{
			    subit = true;
			}
			break;
		    }

		    if (d2 < cosincup2 || cosangle2 < coslimit)
		    {
			if (dislim2)
			{
			    sf_warning ("Hit DIS limit at level %d", pt[count].level[1]);
			}
			else
			{
			    subit = true;
			}
			break;
		    }
		}

		if (subit)
		{
		    if (ncount >= pt_len)
		    {
			morepatches ();
		    }
		    subdivide2 (order, which, cc, &pt[count], &pt[ncount]);
		    ncount++;
		    subdivided = true;
		    if (ncount % 100 == 0)
		    {
			sf_warning ("%d tiles, %f percent", ncount, 100. * (float) count / (float) ncount);
		    }
		}
	    }

	    if (!subdivided)
	    {
		count++;
	    }
	}
	sf_warning ("%d tiles after refinement", ncount);
    }

/*
 * Map over the tiles we used to calculate everything (patch) to
 * the graphical (polypatch) tiles that only contain information
 * needed to do 3D plots.
 */
/* Count how many output tiles are generated */
    npoints = 0;
/*
 * Keep track of how many singular polygons. (They aren't necessarily
 * discarded; it depends on the "singular" option.)
 */
    discarded = 0;
/*
 * Loop over "patch" tiles. These contain info about what's going on
 * in a given Wyoming-shaped piece of Phase space. Each of these tiles
 * generates 3 "polypatch" tiles, one for each wavetype. These may get
 * discarded because of a clipping option, however.
 */
    for (count = 0; count < ncount; count++)
    {
/* Loop over 3 wavetypes */
	for (ii = 0; ii < 3; ii++)
	{
/*
 * Skip surfaces we don't want plotted.
 */
	    if (what[ii] == 0)
		continue;

/*
 * Handle various complicated clipping options for throwing away some
 * surfaces or pieces of surfaces.
 */
	    if (pt[count].clipit && skip != ii && skip < 3)
		continue;

	    if (pt[count].clipit && skip == 3 &&
		!(pt[count].color[ii].red >= pt[count].color[ii].green &&
		  pt[count].color[ii].red >= pt[count].color[ii].blue)
		)
		continue;

	    if (pt[count].clipit && skip == 4 &&
		!(pt[count].color[ii].green >= pt[count].color[ii].red &&
		  pt[count].color[ii].green >= pt[count].color[ii].blue)
		)
		continue;

	    if (pt[count].clipit && skip == 5 &&
		!(pt[count].color[ii].blue >= pt[count].color[ii].green &&
		  pt[count].color[ii].blue >= pt[count].color[ii].red)
		)
		continue;

/*
 * Do group surface
 */
	    if (which)
	    {
/* Find center of tile */
		xnav = 0.;
		ynav = 0.;
		znav = 0.;

		for (jj = 1; jj < 5; jj++)
		{
		    xn = pt[count].cn[jj].im[ii].co[0];
		    yn = pt[count].cn[jj].im[ii].co[1];
		    zn = pt[count].cn[jj].im[ii].co[2];
		    gslo = pt[count].cn[jj].im[ii].gslo;

		    points[npoints].xxn[jj] = xn;
		    points[npoints].yyn[jj] = yn;
		    points[npoints].zzn[jj] = zn;
		    points[npoints].gvel[jj] = 1. / gslo;

		    xnav += .25 * xn / gslo;
		    ynav += .25 * yn / gslo;
		    znav += .25 * zn / gslo;
		}

/*
 * Clip by location of center in Group space. (Last chance to discard
 * a tile because of its color or place of residence.)
 */
		if (or)
		{
		    if ((xnav < xxmin || xnav > xxmax) &&
			(ynav < yymin || ynav > yymax) &&
			(znav < zzmin || znav > zzmax))
			continue;
		}
		else
		{
		    if (xnav < xxmin || xnav > xxmax ||
			ynav < yymin || ynav > yymax ||
			znav < zzmin || znav > zzmax)
			continue;
		}

		points[npoints].xxn[0] = xnav;
		points[npoints].yyn[0] = ynav;
		points[npoints].zzn[0] = znav;

/* Let Vrgraf know this information has been set */
		points[npoints].gvel[0] = true;
	    }
	    else
	    {
		points[npoints].gvel[0] = false;
	    }

/* Save which wavetype this polypatch tile is */
	    points[npoints].kind = ii;

/* Phase surface. Always do this one */
	    {
/* Find center in Phase space... */
		xav = 0.;
		yav = 0.;
		zav = 0.;

		for (jj = 1; jj < 5; jj++)
		{

		    x0 = pt[count].cn[jj].tr[0].co[0];
		    y0 = pt[count].cn[jj].tr[0].co[1];
		    z0 = pt[count].cn[jj].tr[0].co[2];
		    pslo0 = pt[count].cn[jj].tr[0].wv[ii].pslo;

		    points[npoints].xx[jj] = x0;
		    points[npoints].yy[jj] = y0;
		    points[npoints].zz[jj] = z0;
		    points[npoints].pslo[jj] = pslo0;

		    xav += .25 * x0 * pslo0;
		    yav += .25 * y0 * pslo0;
		    zav += .25 * z0 * pslo0;
		}
		points[npoints].xx[0] = xav;
		points[npoints].yy[0] = yav;
		points[npoints].zz[0] = zav;

		points[npoints].pslo[0] = true;
	    }

/*
 * Save Lat and Long of center of tile in Phase space.
 */
	    points[npoints].longi = (360. / TOPI) * atan2 (yav, xav);
	    points[npoints].lati = (360. / TOPI) *
		asin (zav / sqrt (xav * xav + yav * yav + zav * zav));

/*
 * Particle motion direction
 */
	    points[npoints].part[0] = pt[count].cn[0].tr[0].wv[ii].part[0];
	    points[npoints].part[1] = pt[count].cn[0].tr[0].wv[ii].part[1];
	    points[npoints].part[2] = pt[count].cn[0].tr[0].wv[ii].part[2];

/*
 * Whether to draw a given tile border at all (because it was one of
 * the original borders before all the subdividing), or draw it
 * super-fat and White if it happens to have been on an axis plane
 * originally.
 */
	    for (jj = 0; jj < 4; jj++)
	    {
		points[npoints].drawmask[jj] = pt[count].drawmask[jj];
	    }

/* Color (function of particle motion direction) for this wavetype */
	    points[npoints].red = pt[count].color[ii].red;
	    points[npoints].green = pt[count].color[ii].green;
	    points[npoints].blue = pt[count].color[ii].blue;

/*
 * Check to see if this tile is "singular"
 */

/*
 * First sign: particle motion direction at a corner deviates
 * significantly from particle motion direction at center.
 */
	    cosangle = 1.;
	    for (jj = 1; jj < 5; jj++)
	    {
		cosangle1 = fabs (dot (
				      pt[count].cn[0].tr[0].wv[ii].part,
				      pt[count].cn[jj].tr[0].wv[ii].part
				      ));
		if (cosangle1 < cosangle)
		    cosangle = cosangle1;
	    }
/*
 * Second sign: tile is too big in Group space
 */
	    if (which)
	    {
		d1 = dot (pt[count].cn[4].im[ii].co,
			  pt[count].cn[2].im[ii].co);
		d2 = dot (pt[count].cn[1].im[ii].co,
			  pt[count].cn[3].im[ii].co);
	    }
	    else
	    {
/* Fudge it if we didn't calculate Group stuff */
		d1 = d2 = 1.;
	    }

/*
 * If any of those, tag it for possible special treatment, and
 * spew out some info about it to "singout". (It's handy to just
 * run "Vr3d" and have it tell you exactly where all the singularities
 * are in both Phase and Group coordinates!)
 */
	    if (cosangle < coslimit ||
		d1 < cosincup2x2 || d2 < cosincup2x2)
	    {
		points[npoints].singular = true;
		discarded++;

		fprintf (singout, "\nSingular polygon #%d, kind %d, error angle %f\n",
			 discarded,
			 points[npoints].kind,
			 (360. / TOPI) * (float) acos (cosangle));

		sum = 0.;
		if (which)
		{
		    for (jj = 1; jj < 5; jj++)
		    {
			xn = points[npoints].xxn[jj] * points[npoints].gvel[jj] -
			    points[npoints].xxn[0];
			yn = points[npoints].yyn[jj] * points[npoints].gvel[jj] -
			    points[npoints].yyn[0];
			zn = points[npoints].zzn[jj] * points[npoints].gvel[jj] -
			    points[npoints].zzn[0];
			sum += sqrt (xn * xn + yn * yn + zn * zn);
		    }
		}

		if (which)
		{
		    fprintf (singout, "   Lat: %f, Long: %f, Size: %f, Diag: %f, %f\n",
			     points[npoints].lati, points[npoints].longi,
			     (float) (sum / 2.),
			     (360. / TOPI) * (float) acos (d1),
			     (360. / TOPI) * (float) acos (d2));
		}
		else
		{
		    fprintf (singout, "   Lat: %f, Long: %f\n",
			     points[npoints].lati, points[npoints].longi);
		}

		if (which)
		{
		    x = points[npoints].xxn[0];
		    y = points[npoints].yyn[0];
		    z = points[npoints].zzn[0];
		    fprintf (singout, "  Group Coord: x=%f, y=%f, z=%f\n",
			     x, y, z);
		}
/* We can always do this one */
		{
		    x = points[npoints].xx[0];
		    y = points[npoints].yy[0];
		    z = points[npoints].zz[0];
		    fprintf (singout, "  Phase Coord: x=%f, y=%f, z=%f\n",
			     x, y, z);
		}


		for (jj = 1; jj < 5; jj++)
		{
		    if (which)
		    {
			fprintf (singout, "     corner %d:  xr:%f\tyr:%f\tzr:%f\tr:%f\n", jj,
				 (float) points[npoints].xxn[jj],
				 (float) points[npoints].yyn[jj],
				 (float) points[npoints].zzn[jj],
				 (float) points[npoints].gvel[jj]);
		    }
/* We can always do this one */
		    {
			fprintf (singout, "     corner %d:  xw:%f\tyw:%f\tzw:%f\tr:%f\n", jj,
				 (float) points[npoints].xx[jj],
				 (float) points[npoints].yy[jj],
				 (float) points[npoints].zz[jj],
				 (float) points[npoints].pslo[jj]);
		    }

		    fprintf (singout, "      part:     xp:%f\typ:%f\tzp:%f\n",
			     (float) pt[count].cn[jj].tr[0].wv[ii].part[0],
			     (float) pt[count].cn[jj].tr[0].wv[ii].part[1],
			     (float) pt[count].cn[jj].tr[0].wv[ii].part[2]);
		}
	    }
	    else
	    {
		points[npoints].singular = false;
	    }


	    npoints++;
	    if (npoints >= points_len)
	    {
		morepolys ();
	    }
	}
    }
    sf_warning ("%d polygons, %d singular.", npoints, discarded);

    sf_setformat(out,"native_uchar");

    patchlen = sizeof (struct polypatch);
    sf_putint(out,"n1", patchlen);
    sf_putint(out,"n2", npoints);

    raw = sf_ucharalloc(patchlen);

    /* Write it out */
    for (i=0; i < npoints; i++) {
	memcpy(raw, &points[i], patchlen);
	sf_ucharwrite (raw,patchlen, out);
    }

    exit(0);
}

static void morepatches(void)
{
    pt_len += NPOINTS;
    pt = (struct patch *) realloc ((char *) pt, (unsigned) (pt_len * sizeof (struct patch)));
    if (pt == NULL) sf_error("Can't reallocate space for \"pt\"");
}

static void morepolys(void)
{
    points_len += NPATCH;
    points = (struct polypatch *) realloc ((char *) points, (unsigned) (points_len * sizeof (struct polypatch)));
    if (points == NULL) sf_error("Can't reallocate space for \"points\"");
}
