/* 
 * 3-D traveltime interpolation (from rays to Cartesian cube)
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
    int  pick,fill;
    axa az,ax,ay;    /* Cartesian coordinates */
    int iz,ix,iy;
    axa at,ag,ah,aa; /* Ray coordinates */
    int it,ig,ih;
    
    sf_file Fw; /* wavefronfs file */
    sf_file Ft; /* traveltime file */

    float ***tt=NULL; /* traveltime cube */
    float ***ll=NULL; /* ray length cube */

    pt3d   **wm=NULL; /* wavefront it-1 */
    pt3d   **wo=NULL; /* wavefront it   */
    float  **wl=NULL; /* ray length     */

    pt3d Tm,To; /* points on ray ig,ih: it-1,it*/
    vc3d TmTo;

    double l;
    float  t;
    
    
/*------------------------------------------------------------*/

    sf_init(argc,argv);
    if(! sf_getbool("verb",&verb)) verb=false;
    if(! sf_getint ("pick",&pick)) pick=0;
    if(! sf_getint ("fill",&fill)) fill=0;

    /* wavefronts file (a,g,h,t) */
    Fw = sf_input("in");    

    iaxa(Fw,&aa,1); aa.l=" "; if(verb) raxa(aa);
    iaxa(Fw,&ag,2); ag.l="g"; if(verb) raxa(ag);
    iaxa(Fw,&ah,3); ah.l="h"; if(verb) raxa(ah);
    iaxa(Fw,&at,4); at.l="t"; if(verb) raxa(at);

/*------------------------------------------------------------*/

    /* traveltime file (z,x,y) */
    Ft = sf_output("out");

    if(! sf_getint  ("nz",&az.n)) az.n=100;
    if(! sf_getfloat("oz",&az.o)) az.o=0.;
    if(! sf_getfloat("dz",&az.d)) az.d=1.;
    az.l="z";

    if(! sf_getint  ("nx",&ax.n)) ax.n=100;
    if(! sf_getfloat("ox",&ax.o)) ax.o=0.;
    if(! sf_getfloat("dx",&ax.d)) ax.d=1.;
    ax.l="x";

    if(! sf_getint  ("ny",&ay.n)) ay.n=1;
    if(! sf_getfloat("oy",&ay.o)) ay.o=0.;
    if(! sf_getfloat("dy",&ay.d)) ay.d=1.;
    ay.l="y";

    aa.n=1;
    aa.o=0.;
    aa.d=1.;
    aa.l="";
    
    oaxa(Ft,&az,1); if(verb) raxa(az);
    oaxa(Ft,&ax,2); if(verb) raxa(ax);
    oaxa(Ft,&ay,3); if(verb) raxa(ay);
    oaxa(Ft,&aa,4);

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

/*------------------------------------------------------------*/

    /* allocate wavefronts */
    wm =     pt3dalloc2(ag.n,ah.n);
    wo =     pt3dalloc2(ag.n,ah.n);

    /* allocate ray length array */
    wl = sf_floatalloc2(ag.n,ah.n);

    /* initialize wavefronts */
    for( ih=0; ih<ah.n; ih++) {
	for( ig=0; ig<ag.n; ig++) {
	    wm[ih][ig].x=wo[ih][ig].x=0;
	    wm[ih][ig].y=wo[ih][ig].y=0;
	    wm[ih][ig].z=wo[ih][ig].z=0;
	    wm[ih][ig].v=wo[ih][ig].v=0;

	    wl[ih][ig]  =0;
	}
    }

/*------------------------------------------------------------*/

    /* init INT */
    hwt3d_init(az,ax,ay,at,ag,ah);

/*------------------------------------------------------------*/

    it=0;
    pt3dread2(Fw,wm,ag.n,ah.n,3);

    /* LOOP over time */
    for (it=1; it<at.n; it++) {
	if(verb) sf_warning("it=%d",it);

	t = at.o + it * at.d; // traveltime to the current wavefront

	pt3dread2(Fw,wo,ag.n,ah.n,3); // read wavefront it

	for (ih=0; ih<ah.n; ih++) { 
	    for (ig=0; ig<ag.n; ig++) {
		
		Tm = wm[ih][ig];
		To = wo[ih][ig];
		TmTo = vec3d(&Tm,&To);
		wl[ih][ig] += len3d(&TmTo); // update ray length

		l = wl[ih][ig]; // ray length to the current wavefront

		if     (pick==2) hwt3d_lint(tt,ll,To,t,l);
		else if(pick==1) hwt3d_tint(tt,ll,To,t,l);
		else             hwt3d_nint(tt,ll,To,t,l);
	    }
	}

	/* step in time */
	for( ih=0; ih<ah.n; ih++) {
	    for( ig=0; ig<ag.n; ig++) {
		wm[ih][ig] = wo[ih][ig];
	    }
	}

    } // end it 

    /* fill holes */
    if(fill>0) hwt3d_fill(tt,fill);

    /* write traveltime cube */
    sf_floatwrite(tt[0][0],az.n*ax.n*ay.n,Ft);

/*------------------------------------------------------------*/

    free(**tt); free(*tt); free(tt);
    free(**ll); free(*ll); free(ll);
    ;           free(*wm); free(wm);
    ;           free(*wo); free(wo);
    ;           free(*wl); free(wl);
    
    exit (0);
}
