#include <rsf.h>
/*^*/

#include "hwt3d.h"

#include "vecop.h"
/*^*/

static axa az,ax,ay;
static axa at,ag,ah;
static float ***vv;

void hwt3d_init(float*** vv_in    /* velocity */,
		axa      az_in    /* z axis   */,
		axa      ax_in    /* x axis   */,
		axa      ay_in    /* y axis   */,
		axa      at_in    /* t axis   */,
		axa      ag_in    /* g axis   */,
		axa      ah_in    /* h axis   */)
/*< initialize hwt2d >*/
{
    az = az_in;
    ax = ax_in;
    ay = ay_in;

    at = at_in;
    ag = ag_in;
    ah = ah_in;

    vv = vv_in;
}

/*------------------------------------------------------------*/
float hwt3d_getv(pt3d P) 
/*< get velocity from 3-D cube >*/
{
    double  z, x, y;
    double rz,rx,ry;
    int    jz,jx,jy;
    double fz,fx,fy;
    float  v;

    x = SF_MIN(SF_MAX(ax.o,P.x),ax.o+(ax.n-2)*ax.d);
    y = SF_MIN(SF_MAX(ay.o,P.y),ay.o+(ay.n-2)*ay.d);
    z = SF_MIN(SF_MAX(az.o,P.z),az.o+(az.n-2)*az.d);

    rx = (x-ax.o)/ax.d;
    jx = (int)(rx);
    fx =       rx-jx;

    ry = (y-ay.o)/ay.d;
    jy = (int)(ry);
    fy =       ry-jy;

    rz = (z-az.o)/az.d;
    jz = (int)(rz);
    fz =       rz-jz;

    v = vv[ jy  ][ jx  ][ jz  ] * (1-fz)*(1-fx)*(1-fy) +
	vv[ jy  ][ jx  ][ jz+1] * (  fz)*(1-fx)*(1-fy) +
	vv[ jy+1][ jx  ][ jz  ] * (1-fz)*(1-fx)*(  fy) +
	vv[ jy+1][ jx  ][ jz+1] * (  fz)*(1-fx)*(  fy) + 
	vv[ jy  ][ jx+1][ jz  ] * (1-fz)*(  fx)*(1-fy) +
	vv[ jy  ][ jx+1][ jz+1] * (  fz)*(  fx)*(1-fy) +
	vv[ jy+1][ jx+1][ jz  ] * (1-fz)*(  fx)*(  fy) +
	vv[ jy+1][ jx+1][ jz+1] * (  fz)*(  fx)*(  fy);

    return(v);
}

/*------------------------------------------------------------*/

bool hwt3d_cusp(pt3d Tm, /* it-1,ig   ,ih   */
		pt3d To, /* it  ,ig   ,ih   */
		pt3d Gm, /* it  ,ig-1 ,ih   */
		pt3d Gp, /* it  ,ig+1 ,ih   */
		pt3d Hm, /* it  ,ig   ,ih-1 */
		pt3d Hp) /* it  ,ig   ,ih+1 */
/*< find cusp points >*/
{
    int sGmHm, sGmHp;
    int sGpHm, sGpHp;
    int sHmGm, sHmGp;
    int sHpGm, sHpGp;
    bool c;
    
    sGmHm = SF_SIG( jac3d(&To,&Tm,&Gm,&Hm) );
    sGmHp = SF_SIG( jac3d(&To,&Tm,&Gm,&Hp) );

    sGpHm = SF_SIG( jac3d(&To,&Tm,&Gp,&Hm) );
    sGpHp = SF_SIG( jac3d(&To,&Tm,&Gp,&Hp) );

    sHmGm = SF_SIG( jac3d(&To,&Tm,&Hm,&Gm) );
    sHmGp = SF_SIG( jac3d(&To,&Tm,&Hm,&Gp) );

    sHpGm = SF_SIG( jac3d(&To,&Tm,&Hp,&Gm) );
    sHpGp = SF_SIG( jac3d(&To,&Tm,&Hp,&Gp) );

    c = sGmHm*sGmHp==1 ||
	sGpHm*sGpHp==1 ||
	sHmGm*sHmGp==1 ||
	sHpGm*sHpGp==1;

    if(c)
	return true;
    else
	return false;
}

/*------------------------------------------------------------*/

pt3d hwt3d_wfttr(pt3d Tm, 
		 pt3d To,
		 pt3d Gm, 
		 pt3d Gp,
		 pt3d Hm,
		 pt3d Hp)
/*< wavefront tracing >*/
{
    pt3d Tp;
    Tp=hwt3d_step(Tm,To,Gm,Gp,Hm,Hp);
    return(Tp);
}

/*------------------------------------------------------------*/

pt3d hwt3d_raytr(pt3d Tm, 
		 pt3d To)
/*< ray tracing >*/
{
    pt3d Tp;
    pt3d Gm,Gp,Hm,Hp;
    vc3d TmTo;

    double a1,a2,a3;
    vc3d   v1,v2,v3;
    vc3d   qq,uu,ww; /* unit vectors */
    float ro;

/*    printpt3d(Tm);*/
/*    printpt3d(To);*/

    ro = To.v * at.d * 0.1;

    TmTo = vec3d(&Tm,&To);

    v1 = axa3d(1);
    v2 = axa3d(2);
    v3 = axa3d(3);

    a1 = ang3d(&TmTo, &v1); a1 = SF_ABS(a1);
    a2 = ang3d(&TmTo, &v2); a2 = SF_ABS(a2);
    a3 = ang3d(&TmTo, &v3); a3 = SF_ABS(a3);

/*    sf_warning("a1=%g a2=%g a3=%g",a1,a2,a3);*/

/*    if     (a1>=a2 && a1<=a3) ww=v1;*/
/*    else if(a2>=a1 && a2<=a3) ww=v2;*/
/*    else                      ww=v3;*/

    if(      SF_ABS(a1-90) <= SF_ABS(a2-90) &&
	     SF_ABS(a1-90) <= SF_ABS(a3-90) )
	ww=v1;
    else if( SF_ABS(a2-90) <= SF_ABS(a1-90) &&
	     SF_ABS(a2-90) <= SF_ABS(a3-90) )
	ww=v2;
    else 
	ww=v3;

/*    sf_warning("TmTo: %g %g %g",TmTo.dx,TmTo.dy,TmTo.dz);*/
/*    sf_warning(" ww : %g %g %g",  ww.dx,  ww.dy,  ww.dz);*/
/*    sf_warning(" v3 : %g %g %g",  v3.dx,  v3.dy,  v3.dz);*/

    qq   = scl3d(&ww,+1);
    uu   = vcp3d(&TmTo,&qq); /* uu = TmTo x qq */
    qq   = nor3d(&uu);       /* qq = uu / |uu| */
    uu   = scl3d(&qq,ro);    /* uu = qq * ro */
    Gm   = tip3d(&To,&uu);   /* Gm at tip of uu from To */
    Gm.v = hwt3d_getv(Gm);
/*    printpt3d(Gm);*/

    qq   = scl3d(&ww,-1);
    uu   = vcp3d(&TmTo,&qq); 
    qq   = nor3d(&uu);
    uu   = scl3d(&qq,ro);
    Gp   = tip3d(&To,&uu);
    Gp.v = hwt3d_getv(Gp);
/*    printpt3d(Gp);*/

    uu = vec3d(&Gm,&Gp);
    ww = nor3d(&uu);

    qq   = scl3d(&ww,+1);
    uu   = vcp3d(&TmTo,&qq); 
    qq   = nor3d(&uu);
    uu   = scl3d(&qq,ro);
    Hm   = tip3d(&To,&uu);
    Hm.v = hwt3d_getv(Hm);
/*    printpt3d(Hm);*/

    qq   = scl3d(&ww,-1);
    uu   = vcp3d(&TmTo,&qq); 
    qq   = nor3d(&uu);
    uu   = scl3d(&qq,ro);
    Hp   = tip3d(&To,&uu);
    Hp.v = hwt3d_getv(Hp);
/*    printpt3d(Hp);*/

    Tp=hwt3d_step(Tm,To,Gm,Gp,Hm,Hp);
    return(Tp);
}

/*------------------------------------------------------------*/

pt3d hwt3d_step(pt3d Tm, 
		pt3d To, 
		pt3d Gm, 
		pt3d Gp,
		pt3d Hm,
		pt3d Hp)
/*< one HWT time step >*/
{
    double gdx,gdy,gdz,gdr;
    double hdx,hdy,hdz,hdr;
    double ddx,ddy,ddz;
    double ax,bx,ay,by,az,bz;
    double a,b,c,del;

    double dx,dy,dz;

    pt3d Sm,Sp,S;
    vc3d TmTo,ToSm,ToSp;
    double am,ap;
    float ro;

/*------------------------------------------------------------*/

    ro = To.v * at.d;

    gdx = Gp.x-Gm.x;
    gdy = Gp.y-Gm.y;
    gdz = Gp.z-Gm.z;
    gdr =(Gp.v-Gm.v)*at.d;

    hdx = Hp.x-Hm.x;
    hdy = Hp.y-Hm.y;
    hdz = Hp.z-Hm.z;
    hdr =(Hp.v-Hm.v)*at.d;

    ddz = gdy*hdx - gdx*hdy;
    ddx = gdz*hdy - gdy*hdz;
    ddy = gdx*hdz - gdz*hdx;

/*    sf_warning("vvvvvvvvvvvvvvvvvvvvvvvvvv");*/

/*    sf_warning("ddx=%g ddy=%g ddz=%g",ddx,ddy,ddz);*/

    if( SF_ABS(ddz) >=SF_ABS(ddx) && 
	SF_ABS(ddz) > SF_ABS(ddy)) {

/*	sf_warning("use ddz=%g",ddz);*/

	ax  = gdr*hdy - hdr*gdy; 
	bx  = gdy*hdz - gdz*hdy; 
	ay  = gdr*hdx - gdx*hdr; 
	by  = gdx*hdz - gdz*hdx; 
	
	a = ddz*ddz + bx*bx + by*by;
	b =         -(ax*bx + ay*by);           b *= ro;
	c =           ax*ax + ay*ay - ddz*ddz;  c *= ro*ro;
	
	del = SF_MAX(b*b - a*c , 0.);
	del = sqrtf(del);
		
	dz = (-b+del)/a;  Sp.z = To.z + dz;
	dx = ro*ax-dz*bx; Sp.x = To.x + dx/  ddz;
	dy = ro*ay-dz*by; Sp.y = To.y + dy/(-ddz);
	
	dz = (-b-del)/a;  Sm.z = To.z + dz;
	dx = ro*ax-dz*bx; Sm.x = To.x + dx/  ddz;
	dy = ro*ay-dz*by; Sm.y = To.y + dy/(-ddz);
	
    } else if( SF_ABS(ddx) >=SF_ABS(ddy) &&
	       SF_ABS(ddx) > SF_ABS(ddz)) {

/*	sf_warning("use ddx=%g",ddx);*/

	ay  = gdr*hdz - hdr*gdz; 
	by  = gdz*hdx - gdx*hdz; 
	az  = gdr*hdy - gdy*hdr; 
	bz  = gdy*hdx - gdx*hdy; 

	a = ddx*ddx + by*by + bz*bz;
	b =         -(ay*by + az*bz);           b *= ro;
	c =           ay*ay + az*az - ddx*ddx;  c *= ro*ro;
	
	del = SF_MAX(b*b - a*c , 0.);
	del = sqrtf(del);

	dx = (-b+del)/a;  Sp.x = To.x + dx;
	dy = ro*ay-dx*by; Sp.y = To.y + dy/  ddx;
	dz = ro*az-dx*bz; Sp.z = To.z + dz/(-ddx);
	
	dx = (-b-del)/a;  Sm.x = To.x + dx;
	dy = ro*ay-dx*by; Sm.y = To.y + dy/  ddx;
	dz = ro*az-dx*bz; Sm.z = To.z + dz/(-ddx);

    } else {

/*	sf_warning("use ddy=%g",ddy);*/

	az  = gdr*hdx - hdr*gdx; 
	bz  = gdx*hdy - gdy*hdx; 
	ax  = gdr*hdz - gdz*hdr; 
	bx  = gdz*hdy - gdy*hdz; 

	a = ddy*ddy + bz*bz + bx*bx;
	b =         -(az*bz + ax*bx);           b *= ro;
	c =           az*az + ax*ax - ddy*ddy;  c *= ro*ro;
	
	del = SF_MAX(b*b - a*c , 0.);
	del = sqrtf(del);

	dy = (-b+del)/a;  Sp.y = To.y + dy;
	dz = ro*az-dy*bz; Sp.z = To.z + dz/  ddy;
	dx = ro*ax-dy*bx; Sp.x = To.x + dx/(-ddy);
	
	dy = (-b-del)/a;  Sm.y = To.y + dy;
	dz = ro*az-dy*bz; Sm.z = To.z + dz/  ddy;
	dx = ro*ax-dy*bx; Sm.x = To.x + dx/(-ddy);
    }

    TmTo = vec3d( &Tm, &To);
    ToSm = vec3d( &To, &Sm);
    ToSp = vec3d( &To, &Sp);
    
    am = ang3d( &TmTo, &ToSm);
    ap = ang3d( &TmTo, &ToSp);

    if(am<ap) S=Sm;
    else      S=Sp;
    S.v = hwt3d_getv(S);

/*    printpt3d(S);*/
/*    sf_warning("^^^^^^^^^^^^^^^^^^^^^^^^^^");*/

    return S;
}

/*------------------------------------------------------------*/
