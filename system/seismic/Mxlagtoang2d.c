/* SS(x-lag) to angle transformation (PP or PS waves) */

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

/*
 * input:  z-SS(h,z)-x [SS=slant-stack]
 * output: z-a-x
 */

int main(int argc, char* argv[])
{
    bool inv, verb;

    sf_file Fstk=NULL; /*    SS(h,z) file */
    sf_file Fang=NULL; /*     AD CIG file */
    sf_file Fgam=NULL; /* vpvs ratio file */
    sf_file Fdip=NULL; /*  dip field file */

    float **stk=NULL, *gam=NULL, *dip=NULL, **ang=NULL; /* I/O arrays */
    float *tmp=NULL;                     /* mapping arrays */

    sf_axis axz; /* depth axis */
    sf_axis axs; /*    SS axis */
    sf_axis axa; /* angle axis */

    int ext;
    int ncig, icig; /* CIG index */

    int   na;
    float oa,da;
    int   ia,iz,is;
    int   fint;

    fint1 sft;
    float a,t,g,d,n,f;

    /*------------------------------------------------------------*/
    sf_init (argc,argv);

    if (!sf_getbool("verb",&verb)) verb=false;  /* verbosity flag */
    if (!sf_getbool("inv", &inv))   inv=false;  /* inverse transformation flag */

    Fstk = sf_input (  "in"); /* SS(h,z) CIG */
    Fang = sf_output( "out"); /*  AD     CIG */
    Fgam = sf_input ("vpvs"); /*  vpvs ratio */
    Fdip = sf_input ( "dip"); /*  dip  field */

    if (SF_FLOAT != sf_gettype(Fstk)) sf_error("Need float input");

    axz = sf_iaxa(Fstk,1); /* depth axis */
    axs = sf_iaxa(Fstk,2); /*    SS axis */

    ncig = sf_leftsize(Fstk,2); /* number of CIGS to process */

    /* angle axis */
    if (!sf_getint  ("na",&na)) na=sf_n(axs);       
    if (!sf_getfloat("da",&da)) da=1./(sf_n(axs)-1);
    if (!sf_getfloat("oa",&oa)) oa=0.;         
    axa = sf_maxa(na,oa,da);
    sf_oaxa(Fang,axa,2);

    if (!sf_getint("extend",&ext)) ext=4;       /* tmp extension */
    /*------------------------------------------------------------*/

    /* I/O arrays */
    stk = sf_floatalloc2(sf_n(axz),sf_n(axs)); /* SS(h,z) CIG */
    gam = sf_floatalloc (sf_n(axz)         ); /*  vpvs ratio */
    dip = sf_floatalloc (sf_n(axz)         ); /*  dip  field */
    ang = sf_floatalloc2(sf_n(axz),sf_n(axa)); /*      AD CIG */

    /* temp array */
    tmp = sf_floatalloc(sf_n(axs));

    /*------------------------------------------------------------*/
    sft = fint1_init(ext,sf_n(axs),0);

    /*------------------------------------------------------------*/
    for (icig=0; icig < ncig; icig++) { /* loop over CIG */
	if(verb) sf_warning("%d of %d",icig+1,ncig);

	sf_floatread(stk[0],sf_n(axz)*sf_n(axs),Fstk);	
	sf_floatread(gam   ,sf_n(axz)         ,Fgam);
	sf_floatread(dip   ,sf_n(axz)         ,Fdip);

	/*------------------------------------------------------------*/
	for (iz=0; iz < sf_n(axz); iz++) {
	    /* loop over depth */

	    g = gam[iz];
	    d = dip[iz];

	    d*=(g*g-1.);

	    for (is = 0; is < sf_n(axs); is++) { 
		/* loop over slant-stack index */
		tmp[is] = stk[is][iz];
	    }
	    fint1_set(sft,tmp);

	    for (ia=0; ia < sf_n(axa); ia++) {
		a = sf_o(axa)+ia*sf_d(axa);          /* ang */
		t = tanf(a/180*SF_PI);             /* tan */
		
		/*
		 * mapping from tan(a) to slant-stack value (n)
		 */
		n = (4*g*t+d*(t*t+1.)) / ( t*t * (g-1)*(g-1) + (g+1)*(g+1) );

		f = (n - sf_o(axs)) / sf_d(axs);
		fint = f;

		if (fint >= 0 && fint < sf_n(axs)) {
		    ang[ia][iz] = fint1_apply(sft,fint,f-fint,false);
		} else {
		    ang[ia][iz] = 0.;
		}
	    } /* a */

	} /* z */
	/*------------------------------------------------------------*/

	sf_floatwrite(ang[0],sf_n(axz)*sf_n(axa),Fang);
    }

    exit(0);
}
