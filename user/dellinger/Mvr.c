/* Plot impulse responses in 2 dimensions */
/*
 *  This routine getpars for c11, c13, c33, c55, c66, and title
 *   and outputs vplot for a Graph of group velocity and the
 *   dispersion relation for P, Sh, and Sv.
 *  Each of these graphs also have plotted in a dotted line their elliptic
 *   approximation. The graphs are normalized by vertical P elastic constant,
 *   unless norm is specified, in which case they are normalized by that.
 */
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
#include <rsf.h>
#include <rsfplot.h>

#define OFST .7
#define OFST2 1.1
#define NOTHING 10000.

static float vr (float ww, float ss, float dwds);
static float phir (float phiw, float ww, float dwds);
static float partic (int js);

struct polar
{
    float           r;
    float           th;
};

static float ss, phiw;
static float del, eps, zet;
static float degtorad;

int main(int argc, char* argv[])
{
    int             jj;
    float           inc;

    char            *string;
    float           dash, dash2, edash;
    int             color;
    int             fat;

    char           *logfile;
    FILE           *logout;

    float           c11, c13, c33, c44, c66;
    float           alf, bet, gam;
    float           ww, dwds;
    float           normr=0, normd=0, normc, scale;
    float           root, root2=0;
    float           temp;
    struct polar   *vrp, *vrsh, *vrsv, *dp, *dsh, *dsv;
    float          *pp=NULL, *psv=NULL;
    int             num, num2, ii;
    char           *title;
    char            junk[100];
    bool            info;
    bool            disp;
    bool            group;
    bool            ellipse;
    bool            line;
    bool            particle;
    int             part, ipart;
    int             jstart, jend;
    int             decimatep, decimates;
    float           groupscale;
    float           phasescale;

    sf_init(argc,argv);

    if (NULL != (logfile = sf_getstring("log")))
    {
	logout = fopen (logfile, "w");
	if (logout == NULL)
	{
	    sf_warning ("Can't open log file.");
	}
    }
    else
    {	
	logout = (FILE *) NULL;
    }

    if (!sf_getint("color", &color)) color=7;
    if (!sf_getint("fat", &fat)) fat=5;

    if (NULL == (string = sf_getstring("only"))) string="all";
    /* (Pdisp, SVdisp, SHdisp, P, SV, SH) */
    if (!sf_getfloat("dash", &dash)) dash=0.0;
    if (!sf_getfloat("edash", &edash)) edash=0.1;
	/* elliptical approximation dash */

    if (!sf_getfloat("vrscale", &scale)) scale=1.;
    /* scales everything by a factor */

    if (!sf_getfloat("groupscale", &groupscale)) groupscale = 1.;
    /* scales only the group stuff */

    if (!sf_getfloat("phasescale", &phasescale)) phasescale = .5;
    /* scales only the phase stuff */

    if (!sf_getfloat("norm", &normc)) normc = 0.;
    normc = sqrtf(normc);

    inc = (.5);
    if (!sf_getfloat("inc", &inc)) inc=0.5;
    /* increment of phi sub w in degrees */
    if (inc == 0.) sf_error("inc can't be zero!");

    jstart = 0;
    jend = 2;

    if (!sf_getbool("line",&line)) line=false;
    /* if draw lines to indicate some important angles. (Angles at
     which triplication "first" occurs, and angles at which pure P
     and Sv modes exist) */
    if (line) jend = 3;

    if (!sf_getbool("info",&info)) info=true;    
    /* if print in small letters the elastic constants across the top */

    if (!sf_getbool("disp",&disp)) disp=true;    
    /* if n, give phase velocity instead of dispersion relation */

    if (!sf_getbool("group",&group)) group=true;
    /* if n, give group slowness instead of group velocity */

    if (!sf_getbool("particle",&particle)) particle=false;
    /* if show particle motion directions */

    if (!sf_getbool("ellipse", &ellipse)) ellipse=true;
    /* if use elliptic approximation */
    if (!ellipse) jstart = 1;

    if (NULL == (title = sf_getstring("title"))) title=" ";

    num = 360. / fabsf (inc) + 1.;

    if (!sf_getfloat("c11", &c11)) c11=0.;
    if (!sf_getfloat("c13", &c13)) c13=0.;
    if (!sf_getfloat("c33", &c33)) c33=0.;
    if (!sf_getfloat("c55", &c44)) c44=0.;
    if (!sf_getfloat("c66", &c66)) c66=0.;

    if (logout != NULL)
	fprintf (logout, "\tc11=%f  c13=%f  c33=%f\n\tc55=%f  c66=%f\n\n", c11, c13, c33, c44, c66);

    del = c33 - c44;
    eps = c11 - c44;
    zet = c13 + c44;
    alf = del * del;
    bet = -2. * del * (del + eps) + 4. * zet * zet;
    gam = (del + eps) * (del + eps) - 4. * zet * zet;

    ii = -1;

    vrp = (struct polar *) sf_alloc(num, sizeof (struct polar));
    vrsh = (struct polar *) sf_alloc(num, sizeof (struct polar));
    vrsv = (struct polar *) sf_alloc(num, sizeof (struct polar));
    dp = (struct polar *) sf_alloc(num, sizeof (struct polar));
    dsh = (struct polar *) sf_alloc(num, sizeof (struct polar));
    dsv = (struct polar *) sf_alloc(num, sizeof (struct polar));
    if (particle)
    {
	pp  = sf_floatalloc (num);
	psv = sf_floatalloc (num);
    }

    vp_init();
    vp_erase ();

    degtorad = SF_PI / 180.;
    if (strcmp (string, "all") == 0)
    {
	vp_color (7);
	vp_fat (5);
	if (group)
	    vp_text (1750. / 600. - 1.2 - OFST, 5250. / 600. - .6, 7, 90 * (0), "Group Vel.");
	else
	    vp_text (1750. / 600. - 1.2 - OFST, 5250. / 600. - .6, 7, 90 * (0), "Group Slow.");

	if (disp)
	{
	    vp_text (3625. / 600. - 1.5 - OFST, 5250. / 600. - .6, 7, 90 * (0), "Dispersion Rel.");
	}
	else
	{
	    vp_text (3625. / 600. - 1.1 - OFST, 5250. / 600. - .6, 7, 90 * (0), "Phase Vel.");
	}
	if (info)
	{
	    vp_fat (0);

	    sprintf (junk, "c13=%-8.2f", c13);
	    vp_text (500. / 600. - OFST + .4, 1.7 + 5600. / 600., 5, 90 * (0), junk);
	    sprintf (junk, "c11=%-8.2f c33=%-8.2f", c11, c33);
	    vp_text (500. / 600. - OFST + .4, 1.4 + 5600. / 600., 5, 90 * (0), junk);
	    sprintf (junk, "c44=%-8.2f c66=%-8.2f", c44, c66);
	    vp_text (500. / 600. - OFST + .4, 1.1 + 5600. / 600., 5, 90 * (0), junk);

/*
  sprintf (junk, "zet=%-8.2f", zet);
  vp_text (4.+ 500./ 600.- OFST +.4, 1.7 + 5600./ 600., 5,90*( 0), junk);
  sprintf (junk, "eps=%-8.2f", eps);
  vp_text (4.+ 500./ 600.- OFST +.4, 1.4 + 5600./ 600., 5,90*( 0), junk);
  sprintf (junk, "del=%-8.2f", del);
  vp_text (4.+ 500./ 600.- OFST +.4, 1.1 + 5600./ 600., 5,90*( 0), junk);
*/
	}
	vp_fat (7);
	vp_text (500. / 600. - OFST, 5600. / 600., 12, 90 * (0), title);
	vp_fat (5);
	sprintf (junk, "SH");
	vp_text (1. - OFST, 750. / 600., 10, 90 * (0), junk);
	if (del > 0 && eps > 0)
	{
	    sprintf (junk, "SV");
	    vp_text (1. - OFST, 2250. / 600., 10, 90 * (0), junk);
	    vp_text (1. - OFST, 3750. / 600., 10, 90 * (0), "P");
	}
	else
	{
	    vp_text (1. - OFST, 2250. / 600., 10, 90 * (0), "?");
	    vp_text (1. - OFST, 3750. / 600., 10, 90 * (0), "?");
	}

	if (jstart == 0)
	{
	    vp_setdash (&edash, &edash, 1);
	    vp_color (5);
	    vp_fat (3);
	    vp_move (950. / 600., 200. / 600. - OFST2);
	    vp_draw (2100. / 600., 200. / 600. - OFST2);
	    vp_color (7);
	    vp_fat (4);
	    vp_text (2200. / 600., 200. / 600. - OFST2, 5, 90 * (0), "Approximating Ellipse");
	    vp_setdash (&edash, &edash, 0);
	    vp_color (color);
	    vp_fat (fat);
	    vp_move (950. / 600., 450. / 600. - OFST2);
	    vp_draw (2100. / 600., 450. / 600. - OFST2);
	    vp_fat (4);
	    vp_text (2200. / 600., 450. / 600. - OFST2, 5, 90 * (0), "Exact");
	}
    }

    for (jj = jstart; jj < jend; jj++)
    {
	switch (jj)
	{
	    case 0:
		vp_setdash (&edash, &edash, 1);
		vp_color (5);
		vp_fat (3);
		break;
	    case 1:
		if (dash > 0.)
		{
		    dash2 = 2.5 * dash;
		    vp_setdash (&dash, &dash2, 1);
		}
		else
		    vp_setdash (&dash, &dash, 0);
		vp_color (color);
		vp_fat (fat);
		break;
	    case 2:
		vp_setdash (&edash, &edash, 0);
		vp_color (4);
		vp_fat (0);
		break;
	}

	if (jj == 2)
	{
	    num2 = 4;
	}
	else
	    num2 = num;
	for (ii = 0; ii < num2; ii++)
	{
	    if (jj != 2)
	    {
		phiw = ((float) ii + .5) * inc * degtorad;
	    }
	    else
	    {
		switch (ii)
		{
		    case 0:
			if ((del + eps != 0.) && del / (del + eps) >= 0. && del / (del + eps) <= 1.)
			    phiw = asinf (sqrtf (del / (del + eps)));
			else
			    phiw = 0.;
			break;
		    case 1:
			phiw = 0.;
			break;
		    case 2:
			phiw = 90. * degtorad;
			break;
		    case 3:
			if ((del + eps - 2. * zet) != 0. &&
			    (del - zet) / (del + eps - 2. * zet) >= 0. &&
			    (del - zet) / (del + eps - 2. * zet) <= 1.)
			    phiw = asinf (sqrtf ((del - zet) / (del + eps - 2. * zet)));
			else
			    phiw = 0.;
			break;
		}
	    }
	    ss = sinf (phiw);
	    ss = ss * ss;

	    if (jj > 0)
	    {
		/*
		 * Evaluate in such a way as to ensure positivity inside root 
		 */
		temp = (ss * (del + eps) - del);
		root = sqrtf (4. * zet * zet * (1 - ss) * ss + temp * temp);
		if (root == 0.)
		    sf_error ("Square root becomes Zero! (Exact)");
	    }
	    else
	    {
		if (alf == 0.)
		    sf_error ("No elliptic approximation possible!");
		root = sqrtf (alf) + ss * (bet / (2. * sqrtf (alf)));
		root2 = bet / (2. * sqrtf (alf));
	    }

	    /* Sh */
	    if (jj == 1)
	    {
		ww = c44 + (c66 - c44) * ss;
		if (ww <= 0.)
		    sf_error ("Vw for Sh becomes imaginary or infinite! (Exact)");
		dwds = c66 - c44;
		if (group)
		    vrsh[ii].r = vr (ww, ss, dwds);
		else
		    vrsh[ii].r = 1. / vr (ww, ss, dwds);

		vrsh[ii].th = phir (phiw, ww, dwds);
		if (disp)
		{
		    dsh[ii].r = 1. / sqrtf (ww);
		}
		else
		{
		    dsh[ii].r = sqrtf (ww);
		}
		dsh[ii].th = phiw;
	    }

	    /* Sv */
	    ww = (.5) * (c33 + c44 + (c11 - c33) * ss - root);
	    if (ww > 0.)
	    {
		if (jj == 0)
		{
		    dwds = (.5) * ((c11 - c33) - root2);
		}
		else
		{
		    dwds = (.5) *
			((c11 - c33) - .5 * (bet + 2. * ss * gam) / root);
		}
		if (group)
		    vrsv[ii].r = vr (ww, ss, dwds);
		else
		    vrsv[ii].r = 1. / vr (ww, ss, dwds);

		vrsv[ii].th = phir (phiw, ww, dwds);
		if (disp)
		{
		    dsv[ii].r = 1. / sqrtf (ww);
		}
		else
		{
		    dsv[ii].r = sqrtf (ww);
		}
		dsv[ii].th = phiw;
		if (jj > 0 && particle)
		{
		    psv[ii] = partic (-1);
		}
	    }
	    else
	    {
		vrsv[ii].r = 0.;
		dsv[ii].r = 0.;
		vrsv[ii].th = NOTHING;
		dsv[ii].th = NOTHING;
	    }

	    /* P */
	    ww = (.5) * (c33 + c44 + (c11 - c33) * ss + root);
	    if (ww > 0.)
	    {
		if (jj == 0)
		{
		    dwds = (.5) * ((c11 - c33) + root2);
		}
		else
		{
		    dwds = (.5) *
			((c11 - c33) + .5 * (bet + 2. * ss * gam) / root);
		}
		if (group)
		    vrp[ii].r = vr (ww, ss, dwds);
		else
		    vrp[ii].r = 1. / vr (ww, ss, dwds);

		vrp[ii].th = phir (phiw, ww, dwds);
		if (disp)
		{
		    dp[ii].r = 1. / sqrtf (ww);
		}
		else
		{
		    dp[ii].r = sqrtf (ww);
		}
		dp[ii].th = phiw;
		if (jj > 0 && particle)
		{
		    pp[ii] = partic (1);
		}
	    }
	    else
	    {
		vrp[ii].r = 0.;
		dp[ii].r = 0.;
		vrp[ii].th = NOTHING;
		dp[ii].th = NOTHING;
	    }
	}

/* normalize */
/* Passes 1 and 2 had better have same amplitude for phiw=0 ! */
	if (jj < 2)
	{
	    if (normc == 0.)
	    {
		normr = fabsf (vrp[0].r) / scale;
		normd = fabsf (dp[0].r) / scale;
		if (normr == 0. || normd == 0.)
		    sf_error ("Can't normalize!");
	    }
	    else
	    {
		if (group)
		    normr = normc / scale;
		else
		    normr = 1. / (normc * scale);

		if (disp)
		{
		    normd = 1. / (normc * scale);
		}
		else
		{
		    normd = normc / scale;
		}
	    }
	}
	for (ii = 0; ii < num2; ii++)
	{
	    if (jj == 1)
	    {
		vrsh[ii].r /= normr;
		dsh[ii].r /= normd;
	    }
	    if (jj == 2)
	    {
		vrsv[ii].r = 2.2 - .2 * (ii == 3);
		vrp[ii].r = 2.2 - .2 * (ii == 3);
		dsv[ii].r = 4.4 - .4 * (ii == 3);
		dp[ii].r = 4.4 - .4 * (ii == 3);
	    }
	    else
	    {
		vrsv[ii].r /= normr;
		vrp[ii].r /= normr;
		dsv[ii].r /= normd;
		dp[ii].r /= normd;
	    }
	}

	if (logout != NULL)
	{
	    for (ii = 0; ii < num2; ii++)
	    {
		if (jj == 0)
		    fprintf (logout, "P ellipse: %f, %f\t\t%f, %f\n", vrp[ii].th / degtorad, vrp[ii].r, dp[ii].th / degtorad, dp[ii].r);
		else
		    fprintf (logout, "P exact: %f, %f\t\t%f, %f\n", vrp[ii].th / degtorad, vrp[ii].r, dp[ii].th / degtorad, dp[ii].r);
	    }
	    for (ii = 0; ii < num2; ii++)
	    {
		if (jj == 0)
		    fprintf (logout, "SV ellipse: %f, %f\t\t%f, %f\n", vrsv[ii].th / degtorad, vrsv[ii].r, dsv[ii].th / degtorad, dsv[ii].r);
		else
		    fprintf (logout, "SV exact: %f, %f\t\t%f, %f\n", vrsv[ii].th / degtorad, vrsv[ii].r, dsv[ii].th / degtorad, dsv[ii].r);
	    }
	    if (jj == 1)
		for (ii = 0; ii < num2; ii++)
		{
		    fprintf (logout, "SH: %f, %f\t\t%f, %f\n", vrsh[ii].th / degtorad, vrsh[ii].r, dsh[ii].th / degtorad, dsh[ii].r);
		}
	}

	if (jj == 1 && particle)
	    part = 2;
	else
	    part = 1;
	for (ipart = 0; ipart < part; ipart++)
	{
	    if (ipart == 1)
	    {
		if (del > 0 && eps > 0)
		{
		    decimatep = 8;
		    decimates = 20;
		}
		else
		    if (del < 0 && eps < 0)
		    {
			decimatep = 8;
			decimates = 8;
		    }
		    else
		    {
			decimatep = 14;
			decimates = 14;
		    }
		vp_fat (0);
		vp_color (6);
	    }
	    else
	    {
		decimatep = 1;
		decimates = 1;
	    }

	    if (strcmp (string, "all") == 0 || strcmp (string, "SH") == 0)
	    {
		if (jj == 1 && ipart == 0)
		{
		    if (strcmp (string, "all") == 0)
		    {
			vp_orig (1750. / 600. - OFST, 750. / 600.);
		    }
		    else
		    {
			vp_orig (0., 0.);
		    }
		    vp_uorig (0., 0.);
		    vp_scale (groupscale, groupscale);
		    vp_penup ();
		    for (ii = 0; ii < num2; ii++)
		    {
			vp_upendn (vrsh[ii].r * sinf (vrsh[ii].th), vrsh[ii].r * cosf (vrsh[ii].th));
		    }
		}
	    }

	    if (strcmp (string, "all") == 0 || strcmp (string, "SV") == 0)
	    {
		if (strcmp (string, "all") == 0)
		{
		    vp_orig (1750. / 600. - OFST, 2250. / 600.);
		}
		else
		{
		    vp_orig (0., 0.);
		}
		vp_uorig (0., 0.);
		vp_scale (groupscale, groupscale);
		vp_penup ();
		for (ii = 0; ii < num2; ii += decimates)
		{
		    if (jj == 2)
		    {
			vp_penup ();
			vp_upendn (0., 0.);
		    }
		    if (vrsv[ii].th == NOTHING ||
			(vrsv[ii].r > 2.5 && strcmp (string, "all") == 0))
		    {
			vp_penup ();
		    }
		    else
		    {
			if (ipart == 0)
			{
			    vp_upendn (vrsv[ii].r * sinf (vrsv[ii].th), vrsv[ii].r * cosf (vrsv[ii].th));
			}
			else
			{
			    vp_penup ();
			    vp_upendn (vrsv[ii].r * sinf (vrsv[ii].th) + .1 * sinf (psv[ii]), vrsv[ii].r * cosf (vrsv[ii].th) + .1 * cosf (psv[ii]));
			    vp_upendn (vrsv[ii].r * sinf (vrsv[ii].th) - .1 * sinf (psv[ii]), vrsv[ii].r * cosf (vrsv[ii].th) - .1 * cosf (psv[ii]));
			}
		    }
		}
	    }

	    if (strcmp (string, "all") == 0 || strcmp (string, "P") == 0)
	    {
		if (strcmp (string, "all") == 0)
		{
		    vp_orig (1750. / 600. - OFST, 3750. / 600.);
		}
		else
		{
		    vp_orig (0., 0.);
		}
		vp_uorig (0., 0.);
		vp_scale (groupscale, groupscale);
		vp_penup ();
		for (ii = 0; ii < num2; ii += decimatep)
		{
		    if (jj == 2)
		    {
			vp_penup ();
			vp_upendn (0., 0.);
		    }
		    if (vrp[ii].th == NOTHING ||
			(vrp[ii].r > 2.5 && strcmp (string, "all") == 0))
		    {
			vp_penup ();
		    }
		    else
		    {
			if (ipart == 0)
			{
			    vp_upendn (vrp[ii].r * sinf (vrp[ii].th), vrp[ii].r * cosf (vrp[ii].th));
			}
			else
			{
			    vp_penup ();
			    vp_upendn (vrp[ii].r * sinf (vrp[ii].th) + .1 * sinf (pp[ii]), vrp[ii].r * cosf (vrp[ii].th) + .1 * cosf (pp[ii]));
			    vp_upendn (vrp[ii].r * sinf (vrp[ii].th) - .1 * sinf (pp[ii]), vrp[ii].r * cosf (vrp[ii].th) - .1 * cosf (pp[ii]));
			}
		    }
		}
	    }

	    if (strcmp (string, "all") == 0 || strcmp (string, "SHdisp") == 0)
	    {
		if (jj == 1 && ipart == 0)
		{
		    if (strcmp (string, "all") == 0)
		    {
			vp_orig (3625. / 600. - OFST, 750. / 600.);
		    }
		    else
		    {
			vp_orig (0., 0.);
		    }
		    vp_uorig (0., 0.);
		    vp_scale (phasescale, phasescale);
		    vp_penup ();
		    for (ii = 0; ii < num2; ii++)
		    {
			vp_upendn (dsh[ii].r * sinf (dsh[ii].th), dsh[ii].r * cosf (dsh[ii].th));
		    }
		}
	    }

	    if (strcmp (string, "all") == 0 || strcmp (string, "SVdisp") == 0)
	    {
		if (strcmp (string, "all") == 0)
		{
		    vp_orig (3625. / 600. - OFST, 2250. / 600.);
		}
		else
		{
		    vp_orig (0., 0.);
		}
		vp_uorig (0., 0.);
		vp_scale (phasescale, phasescale);
		vp_penup ();
		for (ii = 0; ii < num2; ii += decimates)
		{
		    if (jj == 2)
		    {
			vp_penup ();
			vp_upendn (0., 0.);
		    }
		    if (dsv[ii].th == NOTHING ||
			(dsv[ii].r > 5. && strcmp (string, "all") == 0))
		    {
			vp_penup ();
		    }
		    else
		    {
			if (ipart == 0)
			{
			    vp_upendn (dsv[ii].r * sinf (dsv[ii].th), dsv[ii].r * cosf (dsv[ii].th));
			}
			else
			{
			    vp_penup ();
			    vp_upendn (dsv[ii].r * sinf (dsv[ii].th) + .2 * sinf (psv[ii]), dsv[ii].r * cosf (dsv[ii].th) + .2 * cosf (psv[ii]));
			    vp_upendn (dsv[ii].r * sinf (dsv[ii].th) - .2 * sinf (psv[ii]), dsv[ii].r * cosf (dsv[ii].th) - .2 * cosf (psv[ii]));
			}
		    }
		}
	    }

	    if (strcmp (string, "all") == 0 || strcmp (string, "Pdisp") == 0)
	    {
		if (strcmp (string, "all") == 0)
		{
		    vp_orig (3625. / 600. - OFST, 3750. / 600.);
		}
		else
		{
		    vp_orig (0., 0.);
		}
		vp_uorig (0., 0.);
		vp_scale (phasescale, phasescale);
		vp_penup ();
		for (ii = 0; ii < num2; ii += decimatep)
		{
		    if (jj == 2)
		    {
			vp_penup ();
			vp_upendn (0., 0.);
		    }
		    if (dp[ii].th == NOTHING ||
			(dp[ii].r > 5. && strcmp (string, "all") == 0))
		    {
			vp_penup ();
		    }
		    else
		    {
			if (ipart == 0)
			{
			    vp_upendn (dp[ii].r * sinf (dp[ii].th), dp[ii].r * cosf (dp[ii].th));
			}
			else
			{
			    vp_penup ();
			    vp_upendn (dp[ii].r * sinf (dp[ii].th) + .2 * sinf (pp[ii]), dp[ii].r * cosf (dp[ii].th) + .2 * cosf (pp[ii]));
			    vp_upendn (dp[ii].r * sinf (dp[ii].th) - .2 * sinf (pp[ii]), dp[ii].r * cosf (dp[ii].th) - .2 * cosf (pp[ii]));
			}
		    }
		}
	    }

	}
    }

    exit(0);
}


static float vr (float ww, float ss, float dwds)
{
    float           out;

    out = ww + ss * (1. - ss) * dwds * dwds / ww;
    if (out < 0.)
	sf_error ("Vr becomes imaginary!");
    return sqrtf (out);
}

static float phir (float phiw, float ww, float dwds)
{
    return phiw + atanf (dwds * sinf (phiw) * cosf (phiw) / ww);
}

static float partic (int js)
{
    float           top;
    float           bottom;
    float           angle;

    top = 2. * sqrtf (1 - ss) * sqrtf (ss) * zet;
    bottom = js * sqrtf ((ss * (del + eps) - del) * (ss * (del + eps) - del) - 4 * (ss - 1) * ss * zet * zet)
	+ (1. - ss) * del - eps * ss;
    if (bottom == 0.)
	angle = 90. * degtorad;
    else
	angle = atanf (top / bottom);
    if (sinf (phiw) * cosf (phiw) < 0.)
	angle = -angle;
    return angle;
}
