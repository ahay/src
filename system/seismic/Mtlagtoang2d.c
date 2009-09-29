/*  SS(t-lag) to angle transformation (PP or PS waves) */

/*
  Copyright (C) 2007 Colorado School of Mines

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

int main (int argc, char* argv[])
{
    bool inv, verb;

    sf_file Fstk=NULL; /*    SS(t,z) file */
    sf_file Fang=NULL; /*     AD CIG file */
    sf_file Fgam=NULL; /* vpvs ratio file */
    sf_file Fdip=NULL; /*  dip field file */
    sf_file Fvel=NULL; /*   velocity file */

    float **stk=NULL, *gam=NULL, *dip=NULL, *vel=NULL, **ang=NULL; /* I/O arrays */
    float *tmp=NULL;                           /* mapping arrays */

    sf_axis az; /* depth axis */
    sf_axis as; /*    SS axis */
    sf_axis aa; /* angle axis */

    int ext;
    int ncig, icig; /* CIG index */

    int   na;
    float oa,da;
    int   ia,iz,is;
    int   fint;

    fint1 sft;
    float a,c,g,d,n,f,v,e;

    /*------------------------------------------------------------*/
    sf_init (argc,argv);

    if (!sf_getbool("verb",&verb)) verb=false;  /* verbosity flag */
    if (!sf_getbool("inv", &inv))   inv=false;  /* inverse transformation flag */

    Fstk = sf_input(  "in"); /* SS(h,z) CIG */
    Fang = sf_output("out"); /*  AD     CIG */
    Fgam = sf_input("vpvs"); /*  vpvs ratio */
    Fdip = sf_input( "dip"); /*  dip  field */
    Fvel = sf_input( "vel"); /*   velocity  */

    if (SF_FLOAT != sf_gettype(Fstk)) sf_error("Need float input");

    az=sf_iaxa(Fstk,1);  /* depth axis */
    as=sf_iaxa(Fstk,2);  /*    SS axis */

    ncig = sf_leftsize(Fstk,2); /* number of CIGS to process */

    /* angle axis */
    if (!sf_getint  ("na",&na)) na=sf_n(as);
    if (!sf_getfloat("da",&da)) da=1./(sf_n(as)-1);
    if (!sf_getfloat("oa",&oa)) oa=0.;
    aa = sf_maxa(na,oa,da);
    sf_oaxa(Fang,aa,2);

    if (!sf_getint("extend",&ext)) ext=4;       /* tmp extension */
    /*------------------------------------------------------------*/

    /* I/O arrays */
    stk = sf_floatalloc2(sf_n(az),sf_n(as)); /* SS(t,z) CIG */
    gam = sf_floatalloc (sf_n(az));          /*  vpvs ratio */
    dip = sf_floatalloc (sf_n(az));          /*  dip  field */
    vel = sf_floatalloc (sf_n(az));          /*   velocity  */
    ang = sf_floatalloc2(sf_n(az),sf_n(aa)); /*      AD CIG */

    /* temp array */
    tmp = sf_floatalloc(sf_n(as));

    /*------------------------------------------------------------*/
    sft = fint1_init(ext,sf_n(as),0);

    /*------------------------------------------------------------*/
    for (icig = 0; icig < ncig; icig++) { /* loop over CIG */
	if(verb) sf_warning("%d of %d",icig+1,ncig);

	sf_floatread(stk[0],sf_n(az)*sf_n(as),Fstk);
	sf_floatread(gam   ,sf_n(az)         ,Fgam);
	sf_floatread(dip   ,sf_n(az)         ,Fdip);
	sf_floatread(vel   ,sf_n(az)         ,Fvel);
	
	/*------------------------------------------------------------*/
	for (iz = 0; iz < sf_n(az); iz++) { 
	    /* loop over depth */

	    g = gam[iz];
	    d = dip[iz];
	    v = vel[iz];

	    v = v * sqrtf(1+d*d); /* dip correction */
	    e = (1-g)*(1-g) / (4*g);

	    for (is = 0; is < sf_n(as); is++) { 
		/* loop over slant-stack index */
		tmp[is] = stk[is][iz];
	    }
	    fint1_set(sft,tmp);

	    for (ia=0; ia < sf_n(aa); ia++) {
		a = sf_o(aa)+ia*sf_d(aa);          /* ang */
		c = cosf(a/180*SF_PI);             /* cos */

		/* 
		 * mapping from cos(a) to slant-stack value (n)
		 */
		n = v / sqrtf(g) / sqrtf(c*c + e);

		f = (n - sf_o(as)) / sf_d(as);
		fint = f;

		if (fint >= 0 && fint < sf_n(as)) {
		    ang[ia][iz] = fint1_apply(sft,fint,f-fint,false);
		} else {
		    ang[ia][iz] = 0.;
		}
	    } /* a */

	} /* z */
	/*------------------------------------------------------------*/

	sf_floatwrite(ang[0],sf_n(az)*sf_n(aa),Fang);
    }

    exit (0);
}
