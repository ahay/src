/* Plot impulse responses in 2 dimensions */
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

/*
 * Vrtest.c c11=0 c13=0 c33=0 c55=0 c66=0 m=-1 which=0 francis=n flip=n
 *       inc=.5 norm= groupscale=1. phasescale=1. invert=n
 *
 *
 *  flip --- 
 *  invert --- 
 *  which --- 
 */
/* KEYWORDS: group-velocity phase-velocity anisotropy plot */

#include <rsf.h>
#include <rsfplot.h>

#define DIS (1.e-4)

struct polar
{
    float           r;
    float           th;
};

static struct polar w_group (float phiw);
static float iw_phase_ti (float phiw);
static float iw_phase_francis (float phiw);

static float (*iw_phase)(float);
static float c11, c13, c33, c55;
static int mm;
static float chi, del11, del33;
static float wxnmo, wx, wz, wznmo;

int main(int argc, char* argv[])
{
    bool            francis, invert, which, flip;

    float           normc;
    float           groupscale;
    float           phasescale;
    float           inc, degtorad, phiw;
    float           radius, angle;
    int             num, ii;
    float           xx, yy;
    struct polar    temp;
    float           dash;

    sf_init(argc,argv);

    if (!sf_getfloat("dash", &dash)) dash=0.0;

    if (!sf_getbool("francis",&francis)) francis=false;
    if (!francis)
	iw_phase = iw_phase_ti;
    else
	iw_phase = iw_phase_francis;

    if (!sf_getbool("invert",&invert)) invert=false;
    /* reciprocal of plotting radius */

    if (!sf_getfloat("groupscale", &groupscale)) groupscale = 1.;
    /* scales only the group stuff */

    if (!sf_getfloat("phasescale", &phasescale)) phasescale = 1.;
    /* scales only the phase stuff */

    if (!sf_getfloat("norm", &normc)) normc = 1.;

    inc = (.5);
    if (!sf_getfloat("inc", &inc)) inc=0.5;
    /* increment of phi sub w in degrees */
    if (inc == 0.) sf_error("inc can't be zero!");

    num = (360. / inc);

    if (!sf_getfloat("c11", &c11)) c11=1.;
    if (!sf_getfloat("c13", &c13)) c13=.5;
    if (!sf_getfloat("c33", &c33)) c33=1.;
    if (!sf_getfloat("c55", &c55)) c55=.25;

    if (!sf_getint("m", &mm)) mm=-1;

    if (!sf_getbool("which",&which)) which=false;    
    /* transform from phase to group domain or vice versa */

    c11 /= normc;
    c13 /= normc;
    c33 /= normc;
    c55 /= normc;

    chi = c13 + c55;
    del33 = c33 - c55;
    del11 = c11 - c55;

    if (mm == 1)
    {
	wx = c11;
	wxnmo = c55 + chi * chi / del33;
	wz = c33;
	wznmo = c55 + chi * chi / del11;
    }
    else
    {
	wx = c55;
	wxnmo = c11 - chi * chi / del33;
	wz = c55;
	wznmo = c33 - chi * chi / del11;
    }

    if (!sf_getbool("flip",&flip)) flip=false;    
    /* reciprocal of W's used in Francis' approximation */

    if (flip)
    {
	wx = 1. / wx;
	wxnmo = 1. / wxnmo;
	wz = 1. / wz;
	wznmo = 1. / wznmo;
    }

    vp_init ();
    vp_erase ();

    degtorad = SF_PI / 180.;

    if (dash == 0.)
	vp_setdash (&dash, &dash, 0);
    else
	vp_setdash (&dash, &dash, 1);

    vp_dash (0., 0., 0., 0.);
    vp_color (7);
    vp_fat (1);

    if (!which)
    {
	vp_penup ();

	for (ii = 0; ii < num; ii++)
	{
	    phiw = ii * inc * degtorad;

	    radius = sqrt ((*iw_phase) (phiw));
	    if (invert)
		radius = 1. / radius;
	    angle = phiw;

	    xx = phasescale * radius * sinf (angle);
	    yy = phasescale * radius * cosf (angle);

	    vp_upendn (xx, yy);
	}
    }

    if (which)
    {
	vp_penup ();

	for (ii = 0; ii < num; ii++)
	{
	    phiw = ii * inc * degtorad;

	    temp = w_group (phiw);

	    radius = sqrt (temp.r);
	    if (invert)
		radius = 1. / radius;
	    angle = temp.th;

	    xx = groupscale * radius * sinf (angle);
	    yy = groupscale * radius * cosf (angle);

	    vp_upendn (xx, yy);
	}
    }

    exit(0);
}

static struct polar w_group (float phiw)
{
    float           ww;
    struct polar    out;
    float           dwdphiw;

    dwdphiw = ((1. / (*iw_phase) (phiw + DIS)) - (1. / (*iw_phase) (phiw - DIS))) / (2. * DIS);
    ww = 1. / (*iw_phase) (phiw);

    out.r = ww + dwdphiw * dwdphiw / (4. * ww);
    out.th = phiw + atan (dwdphiw / (2. * ww));

    return out;
}

static float iw_phase_ti (float phiw)
{
    float           ss, cc;
    float           vv;

    ss = sinf (phiw);
    ss = ss * ss;
    cc = cosf (phiw);
    cc = cc * cc;

    vv = .5 *
	(c55 + c11 * ss + c33 * cc + mm *
	 sqrt (
	     (del11 * ss - del33 * cc)
	     *
	     (del11 * ss - del33 * cc)
	     +
	     4. * ss * cc * chi * chi
	     )
	    );

    return (1. / vv);
}

static float iw_phase_francis (float phiw)
{
    float           ss, cc;
    float           vv;

/*
 * Francis measures his angles the other way around!
 */
    ss = cosf (phiw);
    ss = ss * ss;
    cc = sinf (phiw);
    cc = cc * cc;

    vv =
	(wx * wx * wx * cc * cc * cc +
	 wx * wx * (2 * wz + wznmo) * cc * cc * ss +
	 wz * wz * (2 * wx + wxnmo) * cc * ss * ss +
	 wz * wz * wz * ss * ss * ss) /
	((wx * cc + wz * ss) *
	 (wx * cc + wz * ss));

    return (1. / vv);
}
