/* 
 * 3-D ray tracing w/ random shooting directions
 * pcs 2005
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

    sf_axis az,ax,ay; /* Cartesian coordinates */
    int iz,ix,iy;
    sf_axis at;       /* time */
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

    int n;
    float o,d;
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
        if(! sf_getint  ("nt",&n)) n=100;
    if(! sf_getfloat("ot",&o)) o=0;
    if(! sf_getfloat("dt",&d)) d=0.001;
    at=sf_maxa(n,o,d);
    sf_setlabel(at,"t");

    /* velocity file */
    Fv = sf_input ("in");
    az=sf_iaxa(Fv,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az);
    ax=sf_iaxa(Fv,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax);
    ay=sf_iaxa(Fv,3); sf_setlabel(ay,"y"); if(verb) sf_raxa(ay);

    vv=sf_floatalloc3(sf_n(az),sf_n(ax),sf_n(ay)); 
    sf_floatread(vv[0][0],sf_n(az)*sf_n(ax)*sf_n(ay),Fv);

    /* traveltime file */
    Ft = sf_output("out");
    sf_oaxa(Ft,az,1);
    sf_oaxa(Ft,ax,2);
    sf_oaxa(Ft,ay,3);

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

    /* source location */
    if(! sf_getfloat("xsou",&xsou)) xsou=sf_o(ax) + sf_n(ax)*sf_d(ax)/2;
    if(! sf_getfloat("ysou",&ysou)) ysou=sf_o(ay) + sf_n(ay)*sf_d(ay)/2;
    if(! sf_getfloat("zsou",&zsou)) zsou=sf_o(az) + sf_n(az)*sf_d(az)/2;
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
	t = sf_o(at) + it * sf_d(at);

	To.x = Tm.x + Tm.v*sf_d(at) * sin(g)*cos(h);
	To.y = Tm.y + Tm.v*sf_d(at) * sin(g)*sin(h);
	To.z = Tm.z + Tm.v*sf_d(at) * cos(g);
	To.v = hwt3d_getv(vv,To);

	TmTo = vec3d(&Tm,&To);
	l   += len3d(&TmTo);
	if     (pick==2) hwt3d_lint(tt,ll,To,t,l);
	else if(pick==1) hwt3d_tint(tt,ll,To,t,l);
	else             hwt3d_nint(tt,ll,To,t,l);

	for(it=2; it<sf_n(at); it++) {
	    t = sf_o(at) + it * sf_d(at);

	    Tp = hwt3d_raytr(vv,Tm,To,scaleray);
	    Tm = To;
	    To = Tp;

	    TmTo = vec3d(&Tm,&To);
	    l   += len3d(&TmTo); 
	    if     (pick==2) hwt3d_lint(tt,ll,To,t,l);
	    else if(pick==1) hwt3d_tint(tt,ll,To,t,l);
	    else             hwt3d_nint(tt,ll,To,t,l);

	} /* end it */	

    } /* end iray */

    /* fill holes */
    if(fill>0) hwt3d_fill(tt,fill);
    
    /* write traveltime cube */
    sf_floatwrite(tt[0][0],sf_n(az)*sf_n(ax)*sf_n(ay),Ft);
    
/*------------------------------------------------------------*/

    free(**vv); free(*vv); free(vv);
    free(**tt); free(*tt); free(tt);
    free(**ll); free(*ll); free(ll);

    exit (0);
}
