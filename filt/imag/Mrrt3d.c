/* 
 * 3-D ray tracing w/ random shooting directions
 * pcs 2005
 */

/*
  Copyright (C) 2004 University of Texas at Austin
  
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
#include <time.h>

#include <rsf.h>
#include "hwt3d.h"

int main(int argc, char* argv[])
{
    bool verb;
    int  pick; /* method for selecting multi-arrival times */
    int  fill; /* hole-filling parameter */

    int scaleray;
    int nray,jray,iray;

    axa az,ax,ay; /* Cartesian coordinates */
    int iz,ix,iy;
    axa at;       /* time */
    int it;

    float xsou,ysou,zsou; /* source coordinates */

    /* shooting directions */
    float gmin,gmax,g;
    float hmin,hmax,h;

    sf_file Fv; /*   velocity file */
    sf_file Ft; /* traveltime file */

    float ***vv=NULL;   /*   velocity cube */
    float ***tt=NULL,t; /* traveltime cube */
    float ***ll=NULL,l; /* ray length cube */

    pt3d Tm,To,Tp;
    vc3d TmTo;  /* Tm-To vector */

    int  seed;
/*------------------------------------------------------------*/

    sf_init(argc,argv);

    if(! sf_getbool(    "verb",&verb    ))     verb=false;
    if(! sf_getint (    "pick",&pick    ))     pick=2;
    if(! sf_getint (    "fill",&fill    ))     fill=1;
    if(! sf_getint ("scaleray",&scaleray)) scaleray=1.;
    if(! sf_getint (    "nray",&nray    ))     nray=1;
    if(! sf_getint (    "jray",&jray    ))     jray=1;

    if(! sf_getfloat("gmin",&gmin)) gmin=-90;
    if(! sf_getfloat("gmax",&gmax)) gmax=+90;
    if(! sf_getfloat("hmin",&hmin)) hmin=0;
    if(! sf_getfloat("hmax",&hmax)) hmax=180;

    if(verb) sf_warning("gmin=%g gmax=%g",gmin,gmax);
    if(verb) sf_warning("hmin=%g hmax=%g",hmin,hmax);

    /* time axis */
    if(! sf_getint  ("nt",&at.n)) at.n=100;
    if(! sf_getfloat("ot",&at.o)) at.o=0;
    if(! sf_getfloat("dt",&at.d)) at.d=0.001;
    at.l="t"; if(verb) raxa(at);

    /* velocity file */
    Fv = sf_input ("in");
    iaxa(Fv,&az,1); az.l="z"; if(verb) raxa(az);
    iaxa(Fv,&ax,2); ax.l="x"; if(verb) raxa(ax);
    iaxa(Fv,&ay,3); ay.l="y"; if(verb) raxa(ay);

    vv=sf_floatalloc3(az.n,ax.n,ay.n); 
    sf_floatread(vv[0][0],az.n*ax.n*ay.n,Fv);

    /* traveltime file */
    Ft = sf_output("out");
    oaxa(Ft,&az,1);
    oaxa(Ft,&ax,2);
    oaxa(Ft,&ay,3);

    tt=sf_floatalloc3(az.n,ax.n,ay.n); 
    ll=sf_floatalloc3(az.n,ax.n,ay.n); 

    for(iy=0;iy<ay.n;iy++) {
	for(ix=0;ix<ax.n;ix++) {
	    for(iz=0;iz<az.n;iz++) {
		tt[iy][ix][iz] = MISSING;
		ll[iy][ix][iz] = MISSING;
	    }
	}
    }

    /* source location */
    if(! sf_getfloat("xsou",&xsou)) xsou=ax.o + ax.n*ax.d/2;
    if(! sf_getfloat("ysou",&ysou)) ysou=ay.o + ay.n*ay.d/2;
    if(! sf_getfloat("zsou",&zsou)) zsou=az.o + az.n*az.d/2;
    if(verb) sf_warning("xsou=%f ysou=%f zsou=%f",xsou,ysou,zsou);

/*------------------------------------------------------------*/

    /* init random numbers */
    if(! sf_getint("seed",&seed)) seed = time(NULL); /* random seed */
    init_genrand((unsigned long) seed);
    
/*------------------------------------------------------------*/

    /* init HWT */
    hwt3d_rand(az,ax,ay,at);
    
/*------------------------------------------------------------*/
    
    /* LOOP over rays */
    for(iray=0; iray<nray; iray++) {
	if(verb && iray%jray == 0) sf_warning("iray=%d of %d",iray,nray);

	/* init ray */
	g = gmin + genrand_real1() * (gmax-gmin); 
	h = hmin + genrand_real1() * (hmax-hmin);

	g *= SF_PI/180;
	h *= SF_PI/180;
	
	/* trace ray */
	l = 0;
	t = 0;

	Tm.x = xsou;
	Tm.y = ysou;
	Tm.z = zsou;
	Tm.v = hwt3d_getv(vv,Tm);

	if     (pick==2) hwt3d_lint(tt,ll,Tm,t,l);
	else if(pick==1) hwt3d_tint(tt,ll,Tm,t,l);
	else             hwt3d_nint(tt,ll,Tm,t,l);

	it=1;
	t = at.o + it * at.d;

	To.x = Tm.x + Tm.v*at.d * sin(g)*cos(h);
	To.y = Tm.y + Tm.v*at.d * sin(g)*sin(h);
	To.z = Tm.z + Tm.v*at.d * cos(g);
	To.v = hwt3d_getv(vv,To);

	TmTo = vec3d(&Tm,&To);
	l   += len3d(&TmTo);
	if     (pick==2) hwt3d_lint(tt,ll,To,t,l);
	else if(pick==1) hwt3d_tint(tt,ll,To,t,l);
	else             hwt3d_nint(tt,ll,To,t,l);

	for(it=2; it<at.n; it++) {
	    t = at.o + it * at.d;

	    Tp = hwt3d_raytr(vv,Tm,To,scaleray);
	    Tm = To;
	    To = Tp;

	    TmTo = vec3d(&Tm,&To);
	    l   += len3d(&TmTo); 
	    if     (pick==2) hwt3d_lint(tt,ll,To,t,l);
	    else if(pick==1) hwt3d_tint(tt,ll,To,t,l);
	    else             hwt3d_nint(tt,ll,To,t,l);

	} // end it

    } // end iray

    /* fill holes */
    if(fill>0) hwt3d_fill(tt,fill);
    
    /* write traveltime cube */
    sf_floatwrite(tt[0][0],az.n*ax.n*ay.n,Ft);
    
/*------------------------------------------------------------*/

    free(**vv); free(*vv); free(vv);
    free(**tt); free(*tt); free(tt);
    free(**ll); free(*ll); free(ll);

    exit (0);
}
