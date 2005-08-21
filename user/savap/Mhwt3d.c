/* 
 * 3-D Huygens wavefront tracing traveltimes 
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
#include <rsf.h>
#include "hwt3d.h"

int main(int argc, char* argv[])
{
    bool verb;
    axa az,ax,ay;    /* Cartesian coordinates */
    axa at,ag,ah,aa; /* Ray coordinates */
    int it,ig,ih;
    
    float xsou,ysou,zsou; /* source coordinates */

    sf_file Fv; /* velocity file */
    sf_file Fw; /* wavefronfs file */

    float ***vv=NULL; /* velocity       */
    pt3d   **wm=NULL; /* wavefront it-1 */
    pt3d   **wo=NULL; /* wavefront it   */
    pt3d   **wp=NULL; /* wavefront it+1 */

    pt3d       Tm;       /* point  on wft it-1 */
    pt3d Hm,Gm,To,Gp,Hp; /* points on wft it   */
    pt3d       Tp;       /* point  on wft it+1 */
/*------------------------------------------------------------*/

    sf_init(argc,argv);
    if(! sf_getbool("verb",&verb)) verb=false;

    /* velocity file */
    Fv = sf_input ("in");
    iaxa(Fv,&az,1); az.l="z"; if(verb) raxa(az);
    iaxa(Fv,&ax,2); ax.l="x"; if(verb) raxa(ax);
    iaxa(Fv,&ay,3); ay.l="y"; if(verb) raxa(ay);

    vv=sf_floatalloc3(az.n,ax.n,ay.n); 
    sf_floatread(vv[0][0],az.n*ax.n*ay.n,Fv);
    
    /* source location */
    if(! sf_getfloat("xsou",&xsou)) xsou=ax.o + ax.n*ax.d/2;
    if(! sf_getfloat("ysou",&ysou)) ysou=ay.o + ay.n*ay.d/2;
    if(! sf_getfloat("zsou",&zsou)) zsou=az.o + az.n*az.d/2;
    if(verb) fprintf(stderr,"xsou=%f ysou=%f zsou=%f\n",xsou,ysou,zsou);    

    /* time axis */
    if(! sf_getint  ("nt",&at.n)) at.n=100;
    if(! sf_getfloat("ot",&at.o)) at.o=0;
    if(! sf_getfloat("dt",&at.d)) at.d=0.001;
    at.l="t";

    /* shooting angle axis */
    if(! sf_getint  ("ng",&ag.n)) ag.n= 360;
    if(! sf_getfloat("og",&ag.o)) ag.o=-180;
    if(! sf_getfloat("dg",&ag.d)) ag.d= 1;
    ag.l="g";

    if(! sf_getint  ("nh",&ah.n)) ah.n= 360;
    if(! sf_getfloat("oh",&ah.o)) ah.o=-180;
    if(! sf_getfloat("dh",&ah.d)) ah.d= 1;
    ah.l="h";
    if(ay.n==1) { ah.n=1; ah.o=0.; } /* 2-D case */

/*------------------------------------------------------------*/
    aa.n=3;
    aa.o=0.;
    aa.d=1.;
    aa.l="";

    /* wavefronts file (a,g,h,t) */
    Fw = sf_output("out");    
    oaxa(Fw,&aa,1); if(verb) raxa(aa);
    oaxa(Fw,&ag,2); if(verb) raxa(ag);
    oaxa(Fw,&ah,3); if(verb) raxa(ah);
    oaxa(Fw,&at,4); if(verb) raxa(at);

    /* set the output to float */
    sf_putint(Fw,"esize",4);
    sf_settype(Fw,SF_FLOAT);

/*------------------------------------------------------------*/

    /* allocate wavefronts */
    wm = pt3dalloc2(ag.n,ah.n);
    wo = pt3dalloc2(ag.n,ah.n);
    wp = pt3dalloc2(ag.n,ah.n);

    /* initialize wavefronts */
    for( ih=0; ih<ah.n; ih++) {
	for( ig=0; ig<ag.n; ig++) {
	    wm[ih][ig].x=wo[ih][ig].x=wp[ih][ig].x=0;
	    wm[ih][ig].y=wo[ih][ig].y=wp[ih][ig].y=0;
	    wm[ih][ig].z=wo[ih][ig].z=wp[ih][ig].z=0;
	    wm[ih][ig].v=wo[ih][ig].v=wp[ih][ig].v=0;
	}
    }

/*------------------------------------------------------------*/

    /* init HWT */
    hwt3d_init(vv,az,ax,ay,at,ag,ah);

/*------------------------------------------------------------*/

    /* construct it=0 wavefront */
    it=0;
    for( ih=0; ih<ah.n; ih++) {
	for( ig=0; ig<ag.n; ig++) {
	    wm[ih][ig].x=xsou;
	    wm[ih][ig].y=ysou;
	    wm[ih][ig].z=zsou;
	    wm[ih][ig].v=hwt3d_getv(wm[ih][ig]);
	}
    }
    pt3dwrite2(Fw,wm,ag.n,ah.n,3); /* write wavefront it=0 */

/*------------------------------------------------------------*/

    /* construct it=1 wavefront */
    it=1;
    for( ih=0; ih<ah.n; ih++) {
	for( ig=0; ig<ag.n; ig++) {
	    double d,g,h;
	    
	    d = at.d * hwt3d_getv(wm[ih][ig]);
	    g = (ag.o+ig*ag.d) * SF_PI/180;
	    h = (ah.o+ih*ah.d) * SF_PI/180;

	    wo[ih][ig].x=xsou + d*sin(g)*cos(h);
	    wo[ih][ig].y=ysou + d*sin(g)*sin(h);
	    wo[ih][ig].z=zsou + d*cos(g);
	    wo[ih][ig].v=hwt3d_getv(wo[ih][ig]);
	}
    }
    pt3dwrite2(Fw,wo,ag.n,ah.n,3); /* write wavefront it=1 */

/*------------------------------------------------------------*/
    /* LOOP over time */
    for (it=2; it<at.n; it++) {
	if(verb) fprintf(stderr,"it=%d\n",it);

	/* boundaries */
	ig=0;      
	for( ih=0; ih<ah.n; ih++) {
	    wp[ih][ig] = hwt3d_raytr(wm[ih][ig],wo[ih][ig]);
	}
	ig=ag.n-1; 
	for( ih=0; ih<ah.n; ih++) {
	    wp[ih][ig] = hwt3d_raytr(wm[ih][ig],wo[ih][ig]);
	}

	ih=0;      
	for( ig=0; ig<ag.n; ig++) {
	    wp[ih][ig] = hwt3d_raytr(wm[ih][ig],wo[ih][ig]);
	}
	ih=ah.n-1; 
	for( ig=0; ig<ag.n; ig++) {
	    wp[ih][ig] = hwt3d_raytr(wm[ih][ig],wo[ih][ig]);
	}

	for (ih=1; ih<ah.n-1; ih++) { 
	    for (ig=1; ig<ag.n-1; ig++) {

		Tm = wm[ih  ][ig  ];
		To = wo[ih  ][ig  ];  

		Gm = wo[ih  ][ig-1];
		Gp = wo[ih  ][ig+1];

		Hm = wo[ih-1][ig  ];
		Hp = wo[ih+1][ig  ];
		
		if(hwt3d_cusp(Tm,To,Gm,Gp,Hm,Hp)) {
		    Tp = hwt3d_raytr(Tm,To);
		} else {
		    Tp = hwt3d_wfttr(Tm,To,Gm,Gp,Hm,Hp);
		}
		wp[ih][ig] = Tp;
	    }
	}

	/* write wavefront it */
	pt3dwrite2(Fw,wp,ag.n,ah.n,3);
	
	/* step in time */
	for( ih=0; ih<ah.n; ih++) {
	    for( ig=0; ig<ag.n; ig++) {
		wm[ih][ig] = wo[ih][ig];
		wo[ih][ig] = wp[ih][ig];
	    }
	}

    } /* end it */

/*------------------------------------------------------------*/
    exit (0);
}
