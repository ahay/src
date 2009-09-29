/* 
 * 3-D traveltime interpolation (from rays to Cartesian cube)
 * pcs 2005
 */

#include <math.h>
#include <rsf.h>
#include "hwt3d.h"

int main(int argc, char* argv[])
{
    bool verb;
    int  pick; /* method for selecting multi-arrival times */
    int  fill; /* hole-filling parameter */

    sf_axis az,ax,ay;    /* Cartesian coordinates */
    int iz,ix,iy;
    sf_axis at,ag,ah,aa; /* Ray coordinates */
    int it,ig,ih;
    
    sf_file Fw; /* wavefronfs file */
    sf_file Ft; /* traveltime file */

    float ***tt=NULL; /* traveltime cube */
    float ***ll=NULL; /* ray length cube */

    pt3d   **wm=NULL; /* wavefront it-1 */
    pt3d   **wo=NULL; /* wavefront it   */
    float  **wl=NULL; /* ray length     */

    pt3d Tm,To; /* points on ray ig,ih: it-1,it*/
    vc3d TmTo;  /* Tm-To vector */

    double l;
    float  t;

    int n;
    float o,d;


/*------------------------------------------------------------*/

    sf_init(argc,argv);
    if(! sf_getbool("verb",&verb)) verb=false;
    if(! sf_getint ("pick",&pick)) pick=2;
    if(! sf_getint ("fill",&fill)) fill=1;

    /* wavefronts file (a,g,h,t) */
    Fw = sf_input("in");    

    aa=sf_iaxa(Fw,1); sf_setlabel(aa," "); if(verb) sf_raxa(aa);
    ag=sf_iaxa(Fw,2); sf_setlabel(ag,"g"); if(verb) sf_raxa(ag);
    ah=sf_iaxa(Fw,3); sf_setlabel(ah,"h"); if(verb) sf_raxa(ah);
    at=sf_iaxa(Fw,4); sf_setlabel(at,"t"); if(verb) sf_raxa(at);

/*------------------------------------------------------------*/

    /* traveltime file (z,x,y) */
    Ft = sf_output("out");

    if(! sf_getint  ("nz",&n)) n=100;
    if(! sf_getfloat("oz",&o)) o=0.;
    if(! sf_getfloat("dz",&d)) d=1.;
    az=sf_maxa(n,o,d);
    sf_setlabel(az,"z");

    if(! sf_getint  ("nx",&n)) n=100;
    if(! sf_getfloat("ox",&o)) o=0.;
    if(! sf_getfloat("dx",&d)) d=1.;
    ax=sf_maxa(n,o,d);
    sf_setlabel(ax,"x");

    if(! sf_getint  ("ny",&n)) n=1;
    if(! sf_getfloat("oy",&o)) o=0.;
    if(! sf_getfloat("dy",&d)) d=1.;
    ay=sf_maxa(n,o,d);
    sf_setlabel(ay,"y");

    aa=sf_maxa(1,0,1);
    sf_setlabel(aa,"");
    
    sf_oaxa(Ft,az,1); if(verb) sf_raxa(az);
    sf_oaxa(Ft,ax,2); if(verb) sf_raxa(ax);
    sf_oaxa(Ft,ay,3); if(verb) sf_raxa(ay);
    sf_oaxa(Ft,aa,4);

    tt=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay)); 
    ll=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay)); 

    for(iy=0;iy<sf_n(ay);iy++) {
	for(ix=0;ix<sf_n(ax);ix++) {
	    for(iz=0;iz<sf_n(az);iz++) {
		tt[iy][ix][iz] = MISSING;
		ll[iy][ix][iz] = MISSING;
	    }
	}
    }

/*------------------------------------------------------------*/

    /* allocate wavefronts */
    wm =     pt3dalloc2(sf_n(ag),sf_n(ah));
    wo =     pt3dalloc2(sf_n(ag),sf_n(ah));

    /* allocate ray length array */
    wl = sf_floatalloc2(sf_n(ag),sf_n(ah));

    /* initialize wavefronts */
    for( ih=0; ih<sf_n(ah); ih++) {
	for( ig=0; ig<sf_n(ag); ig++) {
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
    t = 0;
    pt3dread2(Fw,wm,sf_n(ag),sf_n(ah),3);

    l = 0;
    for (ih=0; ih<sf_n(ah); ih++) { 
	    for (ig=0; ig<sf_n(ag); ig++) {
		Tm = wm[ih][ig];

		if     (pick==3) hwt3d_xint(tt,ll,Tm,t,l);
		else if(pick==2) hwt3d_lint(tt,ll,Tm,t,l);
		else if(pick==1) hwt3d_tint(tt,ll,Tm,t,l);
		else             hwt3d_nint(tt,ll,Tm,t,l);
	    }
    }

    /* LOOP over time */
    for (it=1; it<sf_n(at); it++) {
	if(verb) sf_warning("it=%d of %d",it+1,sf_n(at));

	t = sf_o(at) + it * sf_d(at); /* traveltime to the current wavefront */
	pt3dread2(Fw,wo,sf_n(ag),sf_n(ah),3); /* read wavefront it */

	for (ih=0; ih<sf_n(ah); ih++) { 
	    for (ig=0; ig<sf_n(ag); ig++) {
		
		Tm = wm[ih][ig];
		To = wo[ih][ig];
		TmTo = vec3d(&Tm,&To);

		wl[ih][ig] += len3d(&TmTo); /* update ray length */
		l = wl[ih][ig]; /* ray length to the current wavefront */

		if     (pick==3) hwt3d_xint(tt,ll,To,t,l);
		else if(pick==2) hwt3d_lint(tt,ll,To,t,l);
		else if(pick==1) hwt3d_tint(tt,ll,To,t,l);
		else             hwt3d_nint(tt,ll,To,t,l);
	    }
	}

	/* step in time */
	for( ih=0; ih<sf_n(ah); ih++) {
	    for( ig=0; ig<sf_n(ag); ig++) {
		wm[ih][ig] = wo[ih][ig];
	    }
	}

    } /* end it */

    /* fill holes */
    if(fill>0) hwt3d_fill(tt,fill);

    /* write traveltime cube */
    sf_floatwrite(tt[0][0],sf_n(az)*sf_n(ax)*sf_n(ay),Ft);

/*------------------------------------------------------------*/

    free(**tt); free(*tt); free(tt);
    free(**ll); free(*ll); free(ll);
    ;           free(*wm); free(wm);
    ;           free(*wo); free(wo);
    ;           free(*wl); free(wl);

    exit (0);
}
