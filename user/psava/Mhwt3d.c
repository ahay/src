/* 
 * 3-D Huygens wavefront tracing traveltimes 
 * pcs 2005
 */

#include <math.h>
#include <rsf.h>
#include "hwt3d.h"

int main(int argc, char* argv[])
{
    bool verb, forceray;
    int        scaleray;
    sf_axis az,ax,ay;    /* Cartesian coordinates */
    sf_axis at,ag,ah,aa; /* Ray coordinates */
    int it,ig,ih;
    
    float xsou,ysou,zsou; /* source coordinates */

    sf_file Fv; /*   velocity file */
    sf_file Fw; /* wavefronfs file */

    float ***vv=NULL; /* velocity       */

    pt3d   **wm=NULL; /* wavefront it-1 */
    pt3d   **wo=NULL; /* wavefront it   */
    pt3d   **wp=NULL; /* wavefront it+1 */

    bool   **kk=NULL; /* cusp flag */

    pt3d Tm,To,Tp; /* points on ray ig,ih:    it-1,it,it+1 */
    pt3d Gm,   Gp; /* points on wavefront it: ig-1,   ig+1 */
    pt3d Hm,   Hp; /* points on wavefront it: ih-1,   ih+1 */
    
    int n;
    float o,d;

/*------------------------------------------------------------*/

    sf_init(argc,argv);
    
    if(! sf_getbool(    "verb",&verb    ))     verb=false;
    if(! sf_getbool("forceray",&forceray)) forceray=false;
    if(! sf_getint ("scaleray",&scaleray)) scaleray=1.;

    /* velocity file */
    Fv = sf_input ("in");
    az=sf_iaxa(Fv,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az);
    ax=sf_iaxa(Fv,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax);
    ay=sf_iaxa(Fv,3); sf_setlabel(ay,"y"); if(verb) sf_raxa(ay);

    vv=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay)); 
    sf_floatread(vv[0][0],sf_n(az)*sf_n(ax)*sf_n(ay),Fv);
    
    /* source location */
    if(! sf_getfloat("xsou",&xsou)) xsou=sf_o(ax) + sf_n(ax)*sf_d(ax)/2;
    if(! sf_getfloat("ysou",&ysou)) ysou=sf_o(ay) + sf_n(ay)*sf_d(ay)/2;
    if(! sf_getfloat("zsou",&zsou)) zsou=sf_o(az) + sf_n(az)*sf_d(az)/2;
    if(verb) sf_warning("xsou=%f ysou=%f zsou=%f",xsou,ysou,zsou);

    /* time axis */
    if(! sf_getint  ("nt",&n)) n=100;
    if(! sf_getfloat("ot",&o)) o=0;
    if(! sf_getfloat("dt",&d)) d=0.001;
    at=sf_maxa(n,o,d);
    sf_setlabel(at,"t");

    /* shooting angle axis */
    if(! sf_getint  ("ng",&n)) n=360;
    if(! sf_getfloat("og",&o)) o=-180;
    if(! sf_getfloat("dg",&d)) d=1;
    ag=sf_maxa(n,o,d);
    sf_setlabel(ag,"g");

    if(! sf_getint  ("nh",&n)) n=360;
    if(! sf_getfloat("oh",&o)) o=-180;
    if(! sf_getfloat("dh",&d)) d=1;
    ah=sf_maxa(n,o,d);
    sf_setlabel(ah,"h");
    if(sf_n(ay)==1) { sf_setn(ah,1); sf_seto(ah,0); } /* 2-D case */

/*------------------------------------------------------------*/
    aa=sf_maxa(3,0,3);
    sf_setlabel(aa,"");

    /* wavefronts file (a,g,h,t) */
    Fw = sf_output("out");    
    sf_oaxa(Fw,aa,1); if(verb) sf_raxa(aa);
    sf_oaxa(Fw,ag,2); if(verb) sf_raxa(ag);
    sf_oaxa(Fw,ah,3); if(verb) sf_raxa(ah);
    sf_oaxa(Fw,at,4); if(verb) sf_raxa(at);

    /* set the output to float */
    sf_putint(Fw,"esize",4);
    sf_settype(Fw,SF_FLOAT);

/*------------------------------------------------------------*/

    /* allocate wavefronts */
    wm = pt3dalloc2(sf_n(ag),sf_n(ah));
    wo = pt3dalloc2(sf_n(ag),sf_n(ah));
    wp = pt3dalloc2(sf_n(ag),sf_n(ah));

    /* allocate cusp flag */
    kk = sf_boolalloc2(sf_n(ag),sf_n(ah));

    /* initialize wavefronts */
    for( ih=0; ih<sf_n(ah); ih++) {
	for( ig=0; ig<sf_n(ag); ig++) {
	    wm[ih][ig].x=wo[ih][ig].x=wp[ih][ig].x=0;
	    wm[ih][ig].y=wo[ih][ig].y=wp[ih][ig].y=0;
	    wm[ih][ig].z=wo[ih][ig].z=wp[ih][ig].z=0;
	    wm[ih][ig].v=wo[ih][ig].v=wp[ih][ig].v=0;

	    kk[ih][ig]  = forceray;
	}
    }

/*------------------------------------------------------------*/

    /* init HWT */
    hwt3d_init(az,ax,ay,at,ag,ah);

/*------------------------------------------------------------*/

    /* construct it=0 wavefront */
    it=0;
    if(verb) sf_warning("it=%d of %d",it+1,sf_n(at));
    for( ih=0; ih<sf_n(ah); ih++) {
	for( ig=0; ig<sf_n(ag); ig++) {
	    wm[ih][ig].x=xsou;
	    wm[ih][ig].y=ysou;
	    wm[ih][ig].z=zsou;
	    wm[ih][ig].v=hwt3d_getv(vv,wm[ih][ig]);
/*	    printpt3d(wm[ih][ig]);*/
	}
    }
    pt3dwrite2(Fw,wm,sf_n(ag),sf_n(ah),3); /* write wavefront it=0 */

/*------------------------------------------------------------*/

    /* construct it=1 wavefront */
    it=1;
    if(verb) sf_warning("it=%d of %d",it+1,sf_n(at));
    for( ih=0; ih<sf_n(ah); ih++) {
	for( ig=0; ig<sf_n(ag); ig++) {
	    double d,g,h;
	    
	    d = sf_d(at) * hwt3d_getv(vv,wm[ih][ig]);
	    g = (sf_o(ag)+ig*sf_d(ag)) * SF_PI/180;
	    h = (sf_o(ah)+ih*sf_d(ah)) * SF_PI/180;

	    wo[ih][ig].x=xsou + d*sin(g)*cos(h);
	    wo[ih][ig].y=ysou + d*sin(g)*sin(h);
	    wo[ih][ig].z=zsou + d*cos(g);
	    wo[ih][ig].v=hwt3d_getv(vv,wo[ih][ig]);
	}
    }
    pt3dwrite2(Fw,wo,sf_n(ag),sf_n(ah),3); /* write wavefront it=1 */

/*------------------------------------------------------------*/
    /* LOOP over time */
    for (it=2; it<sf_n(at); it++) {
	if(verb) sf_warning("it=%d of %d",it+1,sf_n(at));

	/* boundaries */
	ig=0;      
	for( ih=0; ih<sf_n(ah); ih++) {
	    wp[ih][ig] = hwt3d_raytr(vv,wm[ih][ig],wo[ih][ig],scaleray);
	}
	ig=sf_n(ag)-1; 
	for( ih=0; ih<sf_n(ah); ih++) {
	    wp[ih][ig] = hwt3d_raytr(vv,wm[ih][ig],wo[ih][ig],scaleray);
	}

	ih=0;      
	for( ig=0; ig<sf_n(ag); ig++) {
	    wp[ih][ig] = hwt3d_raytr(vv,wm[ih][ig],wo[ih][ig],scaleray);
	}
	ih=sf_n(ah)-1; 
	for( ig=0; ig<sf_n(ag); ig++) {
	    wp[ih][ig] = hwt3d_raytr(vv,wm[ih][ig],wo[ih][ig],scaleray);
	}

	for (ih=1; ih<sf_n(ah)-1; ih++) { 
	    for (ig=1; ig<sf_n(ag)-1; ig++) {
		
		Tm = wm[ih  ][ig  ];
		To = wo[ih  ][ig  ];  
		
		Gm = wo[ih  ][ig-1];
		Gp = wo[ih  ][ig+1];
		
		Hm = wo[ih-1][ig  ];
		Hp = wo[ih+1][ig  ];
		
		/* after a cusp, use normal-wavefront HWT */
		if( kk[ih][ig] == false)
		    kk[ih][ig] = hwt3d_cusp(Tm,To,Gm,Gp,Hm,Hp);
		
		if(kk[ih][ig]) Tp = hwt3d_raytr(vv,Tm,To,scaleray);   /* HWT */
		else           Tp = hwt3d_wfttr(vv,Tm,To,Gm,Gp,Hm,Hp);/* HRT */
		
		wp[ih][ig] = Tp;
	    }
	}
	
	/* write wavefront it */
	pt3dwrite2(Fw,wp,sf_n(ag),sf_n(ah),3);
	
	/* step in time */
	for( ih=0; ih<sf_n(ah); ih++) {
	    for( ig=0; ig<sf_n(ag); ig++) {
		wm[ih][ig] = wo[ih][ig];
		wo[ih][ig] = wp[ih][ig];
	    }
	}

    } /* end it */

    /*------------------------------------------------------------*/

    free(**vv); free(*vv); free(vv);
    ;           free(*wm); free(wm);
    ;           free(*wo); free(wo);
    ;           free(*wp); free(wp);
    ;           free(*kk); free(kk);

    /*------------------------------------------------------------*/


    exit (0);
}
