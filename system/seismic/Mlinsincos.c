/* Solve for angle in equation vx*sin(d) + vy*cos(d) = 1/s0. */
/*
  Copyright (C) 2008 University of Texas at Austin

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
#include "fint1.h"

int main(int argc, char* argv[])
{
    fint1 sft;
    int ext,fint;
    float f;

    int na;                        /* number of output angle values */
    int nvx,nvy;                   /* number of velocity in x, number of velocity in y */
    int npos;                      /* number of points where angles are computed */
    int nt;                        /* number of radial lines or polar angles */
    int nr;                        /* number of radius samples on radial lines */

    float da,oa;                   /* increment in output angle, starting position */
    float dt,ot;                   /* increment in polar angle, starting position */
    float dr,orig;                 /* increment in radius, starting position */
    float dvx,ovx;                 /* increment of velocity in x, starting velocity */
    float dvy,ovy;                 /* increment of velocity in y, starting velocity */
    float s0;                      /* reference slowness */

    float **v;                     /* data on 2D velocity squared grid vx-vy */
    float **r;                     /* inverse polar radius in vx-vy plane */
    float **s;                     /* polar angle in vx-vy plane */
    float **a;                     /* output angle values */
    float *tmp;                    /* inverse radius values on radial line */

    int ivx,ivy;                   /* velocity counters */
    int ia;                        /* output angle counter */
    int it;                        /* polar angle counter */
    int ir;                        /* radius counter */

    float vx,vy;                   /* velocity vx and vy */
    float ri;                      /* inverse radius on radial line */
    float t;                       /* polar angle defining radial line */
    float d;                       /* output angle */

    sf_axis aa,apos;
    sf_file in,out;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getfloat("s0",&s0)) sf_error("Need s0=");
    /* reference slowness */

    if (!sf_getint("na",&na)) sf_error("Need na=");
    /* number of angle values. */
    if (!sf_getfloat("da",&da)) sf_error("Need da=");
    /* angle sampling. */
    if (!sf_getfloat("oa",&oa)) sf_error("Need oa=");
    /* angle origin */

    if (!sf_getint("nt",&nt)) nt=180;
    /* number of polar angle for integration. */
    if (!sf_getfloat("dt",&dt)) dt=2.;
    /* polar angle sampling. */
    if (!sf_getfloat("ot",&ot)) ot=0.;
    /* polar angle origin */

    if (!sf_histint(in,"nvx",&nvx)) sf_error("No nvx= in input");
    if (!sf_histfloat(in,"dvx",&dvx)) sf_error("No dvx= in input");
    if (!sf_histfloat(in,"ovx",&ovx)) sf_error("No ovx= in input");

    if (!sf_getint("nr",&nr)) nr=nvx/2;
    /* number of radius on radial lines */
    if (!sf_getfloat("dr",&dr)) dr=dvx;
    /* radius sampling. */

    orig = 1./s0 + 0.5*dr;
   /* radius greater than or equal to 1/s0 */

    if (!sf_getint("extend",&ext)) ext=4;
    /* tmp extension */

    /* read input file parameters */
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"nvy",&nvy)) sf_error("No nvy= in input");
    if (!sf_histfloat(in,"dvy",&dvy)) sf_error("No dvy= in input");
    if (!sf_histfloat(in,"ovy",&ovy)) sf_error("No ovy= in input");

    if ((nvx != nvy) || (dvx != dvy)) sf_error("Need squared grid for vx-vy plane");

    /* number of points where angles are computed. */
    npos = 1;

    /* output file parameters */ 
    aa = sf_maxa(na,oa,da);
    sf_oaxa(out,aa,1);

    apos = sf_maxa(npos,0,1);
    sf_oaxa(out,apos,2);

    sf_putstring(out,"label1","angle");
    sf_putstring(out,"label2","pos");

    /* memory allocations */
    v = sf_floatalloc2(nvx,nvy);
    r = sf_floatalloc2(nvx,nvy);
    s = sf_floatalloc2(nvx,nvy);
    a = sf_floatalloc2(na,npos);

    tmp = sf_floatalloc(nr);

    /* scaled polar coordinates in 2D velocity grid */
    s[0][0] = 0.;
    r[0][0] = 1000000000.;
    for (ivx = 1; ivx < nvx; ivx++) {
	    vx = ovx + ivx*dvx;
	    r[0][ivx] = 1./(s0*fabsf(vx));
	    s[0][ivx] = atan2f(0,vx);
    }
    for (ivy = 1; ivy < nvy; ivy++) {
	vy = ovy + ivy*dvy;
	for (ivx = 0; ivx < nvx; ivx++) {
	    vx = ovx + ivx*dvx;
	    r[ivy][ivx] = 1./(s0*sqrtf(vx*vx + vy*vy));
	    s[ivy][ivx] = atan2f(vy,vx);
	}
    }

    /* read data in velocity array */
    sf_floatread(v[0],nvx*nvy,in);

    /* initialize */
    for (ia = 0; ia < na; ia++) {
	a[0][ia] = 0.;
    }
    sft = fint1_init(ext,nr, 0);

    /* inverse radius values on radial line */
    for (ir = 0; ir < nr; ir++) {
	tmp[ir] = 1./(orig + ir*dr);
    }

    /* Loop on polar angle directions */
    for (it = 0; it < nt; it++) {
	t = ot + it*dt;

	for (ia = 0; ia < na; ia++) {
	    d = oa + ia*da;

            /* cos(d-t) for radial line equals inverse radius */
	    ri = cosf((d-t)/180.*SF_PI);

	    f = (ri - 1/orig)*dr;
	    fint = f;

            /* cumulative sum */
	    if (fint >= 0 && fint < nr) {
		a[0][ia] += fint1_apply(sft,fint,f-fint,false);
	    }
	}
    }
	
    /* output on angle grid */
    sf_floatwrite(a[0],na*npos,out);

    exit(0);
}
