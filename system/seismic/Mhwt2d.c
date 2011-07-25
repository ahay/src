/* 2-D Huygens wavefront tracing traveltimes */

/*
  Copyright (C) 2006 Colorado School of Mines
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
#include <rsf.h>
#include "hwt2d.h"

int main (int argc, char *argv[])
{
    bool verb;
    bool rays;

    sf_axis az,ax;  /* Cartesian coordinates */
    sf_axis at,ag;  /* Ray coordinates */
    int it,ig;
    int nz,nx,nt,ng;
    float dt,dg,ot,og;

    float xsou,zsou; /* source coordinates */

    sf_file Fv=NULL; /* velocity file */
    sf_file Fw=NULL; /* wavefronfs file */

    float **vv=NULL; /* velocity       */
    pt2d   *wm=NULL; /* wavefront it-1 */
    pt2d   *wo=NULL; /* wavefront it   */
    pt2d   *wp=NULL; /* wavefront it+1 */

    pt2d    Ro;    /* point  on wft it-1 */
    pt2d Pm,Po,Pp; /* points on wft it   */
    pt2d    Qo;    /* point  on wft it+1 */

    /*------------------------------------------------------------*/
    sf_init(argc,argv);
    if(! sf_getbool("verb",&verb)) verb=false;
    if(! sf_getbool("rays",&rays)) rays=false;

    /* velocity file */
    Fv = sf_input ("in");
    az = sf_iaxa(Fv,1); sf_setlabel(az,"z"); nz=sf_n(az); if(verb) sf_raxa(az);
    ax = sf_iaxa(Fv,2); sf_setlabel(ax,"x"); nx=sf_n(ax); if(verb) sf_raxa(ax);

    vv=sf_floatalloc2(nz,nx); 
    sf_floatread(vv[0],nz*nx,Fv);

    /* source location */
    if(! sf_getfloat("xsou",&xsou)) xsou=sf_o(ax) + nx*sf_d(ax)/2;
    if(! sf_getfloat("zsou",&zsou)) zsou=sf_o(az) + nz*sf_d(az)/2;
    if(verb) fprintf(stderr,"xsou=%f zsou=%f\n",xsou,zsou);

    /* time axis */
    if(! sf_getint  ("nt",&nt)) nt=100;
    if(! sf_getfloat("ot",&ot)) ot=0;
    if(! sf_getfloat("dt",&dt)) dt=0.001;
    at = sf_maxa(nt,ot,dt); sf_setlabel(at,"t");

    /* shooting angle axis */
    if(! sf_getint  ("ng",&ng)) ng= 360;
    if(! sf_getfloat("og",&og)) og=-180;
    if(! sf_getfloat("dg",&dg)) dg= 1;
    ag = sf_maxa(ng,og,dg); sf_setlabel(ag,"g");

    /*------------------------------------------------------------*/
    /* wavefronts file (g,t) */
    Fw = sf_output("out");
    sf_oaxa(Fw,ag,1); if(verb) sf_raxa(ag);
    sf_oaxa(Fw,at,2); if(verb) sf_raxa(at);

    /* set the output to complex */
    sf_putint(Fw,"esize",8);
    sf_settype(Fw,SF_COMPLEX);

    /*------------------------------------------------------------*/
    /* allocate wavefronts */
    wm = pt2dalloc1(ng);
    wo = pt2dalloc1(ng);
    wp = pt2dalloc1(ng);

    /* initialize wavefronts */
    for( ig=0; ig<ng; ig++) {
	wm[ig].x=wo[ig].x=wp[ig].x=0;
	wm[ig].z=wo[ig].z=wp[ig].z=0;
	wm[ig].v=wo[ig].v=wp[ig].v=0;
    }

    /*------------------------------------------------------------*/
    /* init HWT */
    hwt2d_init(vv,az,ax,at,ag);

    /*------------------------------------------------------------*/
    /* construct it=0 wavefront */
    it=0;
    for( ig=0; ig<ng; ig++) {
	wm[ig].x=xsou;
	wm[ig].z=zsou;
	wm[ig].v=hwt2d_getv(wm[ig]);
    }
    pt2dwrite1(Fw,wm,ng,2); /* write wavefront it=0 */

    /*------------------------------------------------------------*/
    /* construct it=1 wavefront */
    it=1;
    for( ig=0; ig<ng; ig++) {
	double d,g;

	d = dt * hwt2d_getv(wm[ig]);
	g = (og+ig*dg) * SF_PI/180;

	wo[ig].x=xsou + d*sin(g);
	wo[ig].z=zsou + d*cos(g);
	wo[ig].v=hwt2d_getv(wo[ig]);
    }
    pt2dwrite1(Fw,wo,ng,2); /* write wavefront it=1 */

    /*------------------------------------------------------------*/
    /* LOOP over time */
    for (it=2; it<nt; it++) {
	if(verb) fprintf(stderr,"it=%d\n",it);
	
	if(ng>3 && !rays) {
	    /* boundaries */
	    ig=0;      wp[ig] = hwt2d_raytr(wm[ig],wo[ig]);
	    ig=ng-1; wp[ig] = hwt2d_raytr(wm[ig],wo[ig]);

	    for (ig=1; ig<ng-1; ig++) {
		
		Pm = wo[ig-1];
		Po = wo[ig  ];  Qo = wm[ig];
		Pp = wo[ig+1];
		
		if(hwt2d_cusp(Qo,Pm,Po,Pp)) {
		    Ro = hwt2d_raytr(Qo,   Po   );
		} else {
		    Ro = hwt2d_wfttr(Qo,Pm,Po,Pp);
		}
		wp[ig] = Ro;
	    }
	} else {
	    for (ig=0; ig<ng; ig++) {
		Po = wo[ig];  
		Qo = wm[ig];
		Ro = hwt2d_raytr(Qo,Po);
		wp[ig] = Ro;
	    }
	}

	/* write wavefront it */
	pt2dwrite1(Fw,wp,ng,2);

	/* step in time */
	for( ig=0; ig<ng; ig++) {
	    wm[ig] = wo[ig];
	    wo[ig] = wp[ig];
	}
    } /* end it */

    /*------------------------------------------------------------*/

    exit (0);
}
