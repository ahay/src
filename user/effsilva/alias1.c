#include <math.h>
#include <rsf.h>

#include "hwt2d.h"

#include "alias1.h"

/****************************************************************************/
#ifndef _alias1_h

typedef struct boxraysstruct boxrays;
/*^*/

struct boxraysstruct{ 
    float og;  float dg; int ng;
    float ot;  float dt; int nt; 
};
/*^*/

#endif

/*-----------------------------------------------------------------------------*/
int indexMap(
    sf_complex **rays,
    boxrays b, 
    float **vel, 
    sf_axis ax, 
    sf_axis az, 
    float dx, 
    float dz,
    int nf, 
    float df, 
    float fmax, 
    int **map )
/*< index >*/
{
    int flag=0, it, ig, norder=6, index;
    float *x, *z, *dxdg, *dzdg, zero=.02, tg, tg90=9999.;
    float fx, fz, fc, arg, vv;
    pt2d p;
    
    x    = (float *) malloc(b.ng*sizeof(float));
    z    = (float *) malloc(b.ng*sizeof(float));
    dxdg = (float *) malloc(b.ng*sizeof(float));
    dzdg = (float *) malloc(b.ng*sizeof(float));
    
    sf_deriv_init(b.ng, norder, 0.);
    hwt2d_init(vel,az,ax,az,ax);
    
    for(it=0;it<b.nt;it++){
	for(ig=0;ig<b.ng;ig++){
	    x[ig]=crealf(rays[it][ig]);
	    z[ig]=cimagf(rays[it][ig]);
	}
	sf_deriv(x,dxdg);
	sf_deriv(z,dzdg);
	for(ig=0;ig<b.ng;ig++){
	    if(dxdg[ig]>zero) {tg=fabs(dzdg[ig]/dxdg[ig]);
	    }else{tg=tg90;}
	    arg = sqrt(1+tg*tg);
	    p.x = crealf(rays[it][ig]);
	    p.z = cimagf(rays[it][ig]);
	    vv  = hwt2d_getv(p);
	    fx  = vv/(2*dx*arg);
	    /* REVER fz
	       if(tg==tg90) fz=vv/(2*dz);
	       else fz=vv*tg/(2*dz*arg);*/
	    fz=vv/(2*dz);
	    if(fx<fz) fc=fx;
	    else fc=fz;
	    /* DUVIDA
	       if(fc>(fmax-df)) index=0;
	       else index = floor((fmax-df-fc)/df);*/
	    if(fc>fmax) index=0;
	    else index = floor((fmax-fc)/df);
	    if(index<0) index=0;
	    if(index>(nf-1)) index=nf-1;
	    map[it][ig]=index;
	}
    }
    
    sf_deriv_close();
    free(x); free(z); free(dxdg); free(dzdg);
    flag = 1;
    
    return flag;
}

/*-----------------------------------------------------------------------------*/
int angleMap(
    sf_complex **rays,  
    boxrays b, 
    float **map )
/*< angle >*/
{
    int flag=0, it, ig, norder=6;
    float *x, *z, *dxdg, *dzdg, zero=.00001, ang, halfpi;
    
    x    = (float *) malloc(b.ng*sizeof(float));
    z    = (float *) malloc(b.ng*sizeof(float));
    dxdg = (float *) malloc(b.ng*sizeof(float));
    dzdg = (float *) malloc(b.ng*sizeof(float));
    
    sf_deriv_init(b.ng, norder, 0.);
    
    halfpi = 2*atan(1);
    for(it=0;it<b.nt;it++){
	for(ig=0;ig<b.ng;ig++){
	    x[ig]=crealf(rays[it][ig]);
	    z[ig]=cimagf(rays[it][ig]);
	}
	sf_deriv(x,dxdg);
	sf_deriv(z,dzdg);    
	for(ig=0;ig<b.ng;ig++){
	    if(dxdg[ig]>zero) {ang=atan(dzdg[ig]/dxdg[ig]);
	    }else{ang=halfpi;}
	    map[it][ig]=fabs(ang);
	}
    }
    
    sf_deriv_close();
    free(x); free(z); free(dxdg); free(dzdg);
    flag = 1;
    
    return flag;
}

/*-----------------------------------------------------------------------------*/
float isocAngle( 
    float x, 
    float z, 
    float v, 
    float h, 
    float tau)
/*< isoc >*/
{
    float h2, a, a2, ang, halfpi=1.5707963;
    
    if(z>0) {
	h2=h*h; a=v*tau/2; a2=a*a; 
	ang = atan((1-h2/a2)*x/z);
    } else { 
	ang=halfpi;
    }
    
    return ang;
}

/*-----------------------------------------------------------------------------*/
int isocIndex( 
    float x,  float z,  float v, float h,  float tau,
    float dx, float dz, int nf,  float df, float fmax )
/*< isoc index >*/
{
    int index;
    float h2, a, a2, tg, tg90=9999., arg, fx, fz, fc;
    
    if(z>0) {
	h2=h*h; a=v*tau/2; a2=a*a;
	tg = (1-h2/a2)*x/z;
    }else{ tg=tg90;}
    arg = sqrt(1+tg*tg);
    fx  = v/(2*dx*arg);
    fz  = v/(2*dz);
    if(fx<fz) fc=fx;
    else fc=fz;
    /* DUVIDA
       if(fc>(fmax-df)) index=0;
       else index = floor((fmax-df-fc)/df);*/
    if(fc>fmax) index=0;
    else index = floor((fmax-fc)/df);
    
    if(index<0) index=0;
    if(index>(nf-1)) index=nf-1;
    
    return index;
}

/*-----------------------------------------------------------------------------*/
int indexMap2( 
    sf_complex **rays,  
    boxrays b, 
    float **vel, 
    sf_axis ax, 
    sf_axis az, 
    float dx, 
    float dz,
    int nf, 
    float df, 
    float fmax, 
    int **map )
/*< index 2 >*/
{
/* this functions was not tested */
    int flag=0, it, ig, norder=6, index;
    float *x, *z, *dxdt, *dzdt, zero=.02, tg, tg90=9999.;
    float fx, fz, fc, arg, vv;
    pt2d p;
    
    x    = (float *) malloc(b.nt*sizeof(float));
    z    = (float *) malloc(b.nt*sizeof(float));
    dxdt = (float *) malloc(b.nt*sizeof(float));
    dzdt = (float *) malloc(b.nt*sizeof(float));
    
    sf_deriv_init(b.nt, norder, 0.);
    hwt2d_init(vel,az,ax,az,ax);
    
    for(ig=0;ig<b.ng;ig++){
	for(it=0;it<b.nt;it++){
	    x[it]=crealf(rays[it][ig]);
	    z[it]=cimagf(rays[it][ig]);
	}
	sf_deriv(x,dxdt);
	sf_deriv(z,dzdt);
	for(it=0;it<b.nt;it++){
	    if(dxdt[it]>zero) {tg=fabs(dzdt[it]/dxdt[it]);
	    }else{tg=tg90;}
	    arg = sqrt(1+tg*tg);
	    p.x = x[it];
	    p.z = z[it];
	    vv  = hwt2d_getv(p);
	    fx  = vv/(2*dx*arg);
	    /* REVER fz
	       if(tg>zero) fz=tg*vv/(2*dz*arg);
	       else fz=vv/(2*dz);*/
	    fz=vv*tg/(2*dz*arg);
	    if(fx<fz) fc=fx;
	    else fc=fz;
	    if(fz<fx)  fprintf(stderr,"theta=%f fx=%f fz=%f fc=%f\n",57.2957795*atan(tg),fx,fz,fc);
	    if(fc>(fmax-df)) index=0;
	    else index = floor((fmax-df-fc)/df);
	    map[it][ig]=index;
	}
    }
    
    sf_deriv_close();
    free(x); free(z); free(dxdt); free(dzdt);
    flag = 1;
    
    return flag;
}
