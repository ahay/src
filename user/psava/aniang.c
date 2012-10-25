/*
  Copyright (C) 2010 Colorado School of Mines

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
#include <rsf.h>

#include "aniang.h"


/*------------------------------------------------------------*/
float* psitti(vc3d *nn,
	      vc3d *qq,
	      vc3d *tt,
	      vc3d *aa,
	      float tht,   /* reflection angle */
	      float phi,   /* reflection azimuth */
	      sf_axis aps, /* psi axis */
	      float vep,   /* P velocity */
	      float ves,   /* S velocity */
	      float eps,   /* epsilon */
	      float dlt)   /* delta   */
/*< deviation angle in TTI media >*/
{
    float *ttipar;
/*    vc3d bb; */
    int ips;
    float psi,v_s=0.,v_r=0.; /* output variables */
    float ps_,vs_,vr_; /*  local variables */
    float     ts_,tr_;
/*    float     t_s,t_r; */
    float     as_,ar_;
/*    float     a_s,a_r; */
    vc3d  ns_,nr_;
    float snell,snell_;

/*
    float deg2rad,rad2deg;
*/

/*    sf_warning("tht=%g",tht/SF_PI*180.);*/
/*    sf_warning("eps=%g",eps);*/
/*    sf_warning("dlt=%g",dlt);*/
/*    sf_warning("vep=%g",vep);*/
/*    sf_warning("ves=%g",ves);*/

/*
    deg2rad = SF_PI/180;
    rad2deg = 180/SF_PI;
*/

    ttipar=sf_floatalloc(3);
/*    bb=vcp3d(qq,nn); */
    
    snell=1e6;
    psi=0;
    for(ips=0; ips<sf_n(aps); ips++) {
	ps_=(sf_o(aps)+ips*sf_d(aps))/180.*SF_PI;

/*	ts_ = tht-ps_;*/
/*	tr_ = tht+ps_;*/
/*	ns_ = rot3d(&bb,nn,-SF_ABS(ts_)); */
/*	nr_ = rot3d(&bb,nn,+SF_ABS(tr_));*/

	ts_ = ps_-tht;
	tr_ = ps_+tht;
/*	sf_warning("%f %f %f",ps_*180/SF_PI,ts_*180/SF_PI,tr_*180/SF_PI);*/

	ns_ = polvect(phi,ts_,nn,aa);
	nr_ = polvect(phi,tr_,nn,aa);
/*	vout("ns",&ns_);*/
/*	vout("nr",&nr_);*/

	as_ = ang3d(&ns_,tt);
	ar_ = ang3d(&nr_,tt);

	vs_= anivel(as_,vep,ves,eps,dlt);
	vr_= anivel(ar_,vep,ves,eps,dlt);

	snell_ = SF_ABS( SF_ABS(vr_*sin(ts_)) - SF_ABS(vs_*sin(tr_)) );

	if(snell_ < snell) {
	    snell=snell_;
	    psi=ps_;
	    
/*
	    t_s=ts_;
	    t_r=tr_;

	    a_s=as_;
	    a_r=ar_;
*/
	    v_s=vs_;
	    v_r=vr_;
	}
    }

    ttipar[0]=psi/SF_PI*180;
    ttipar[1]=v_s;
    ttipar[2]=v_r;

/*    vout("n",nn);*/
/*    vout("t",tt);*/
/*    vout("q",qq);*/
/*    vout("b",&bb);*/
    
/* 
    sf_warning("ph=%4.1f th=%4.1f ps=%4.1f (ts=%4.1f tr=%4.1f) (as=%4.1f ar=%4.1f) (vs=%5.3f vr=%5.3f)",
	       phi/SF_PI*180,tht/SF_PI*180,psi/SF_PI*180,
	       t_s/SF_PI*180,t_r/SF_PI*180,
	       a_s,a_r,
	       v_s,v_r);
*/

    return ttipar;
}

/*------------------------------------------------------------*/
void vout(char *c,
	  vc3d *v)
{
    sf_warning("%s={%8.4g,%8.4g,%8.4g}",c,v->dx,v->dy,v->dz);
}


/*------------------------------------------------------------*/
float anivel(float alp, /* alpha */
	     float vep, /* P velocity */
	     float ves, /* S velocity */
	     float eps, /* epsilon */
	     float dlt) /* delta */
/*< return TTI velocity as a function of angle >*/
{
    float vel,f;
    float sia,coa;
    float tmp;

    f = 1 - (ves*ves) / (vep*vep);

    sia = sinf(alp/180.*SF_PI); sia=sia*sia;
    coa = cosf(alp/180.*SF_PI); coa=coa*coa;

    tmp = 1 +
	(4*sia/f)*( 2*dlt*coa - eps*(coa-sia) ) + 
	4*eps*eps*sia*sia/(f*f);

    tmp = 1+eps*sia - 0.5*f *(1-sqrt(tmp));

    vel = vep * sqrt(tmp);
    return vel;
}


/*------------------------------------------------------------*/
vc3d rot3dfix(vc3d *nn,
	      vc3d *aa,
	      float phi)
{
    vc3d qq;
    float Qxx,Qxy,Qxz, Qyx,Qyy,Qyz, Qzx,Qzy,Qzz;
    float    nx,ny,nz;
    float    ax,ay,az;
    float    qx,qy,qz;
    float sif,cof;

    nx=nn->dx; ny=nn->dy; nz=nn->dz;
    ax=aa->dx; ay=aa->dy; az=aa->dz;

    sif=sinf(phi);
    cof=cosf(phi);

    Qxx=nx*nx+(ny*ny+nz*nz)*cof;
    Qxy=nx*ny*(1-cof)-nz*sif;
    Qxz=nx*nz*(1-cof)+ny*sif;
    
    Qyx=ny*nx*(1-cof)+nz*sif;
    Qyy=ny*ny+(nz*nz+nx*nx)*cof;
    Qyz=ny*nz*(1-cof)-nx*sif;
    
    Qzx=nz*nx*(1-cof)-ny*sif;
    Qzy=nz*ny*(1-cof)+nx*sif;
    Qzz=nz*nz+(nx*nx+ny*ny)*cof;

    qx = Qxx*ax + Qxy*ay + Qxz*az;
    qy = Qyx*ax + Qyy*ay + Qyz*az;
    qz = Qzx*ax + Qzy*ay + Qzz*az;

    qq.dx=qx;
    qq.dy=qy;
    qq.dz=qz;
    
    return qq;
}

/*------------------------------------------------------------*/
vc3d polvect(float phi,
	     float tht,
	     vc3d *nn,
	     vc3d *aa)
/*<  >*/
{
    vc3d jj;
    vc3d nj;

    jj = rot3d(nn,aa,phi+0.5*SF_PI);

/*    vout("nn",nn);*/
/*    vout("aa",aa);*/
/*    sf_warning("phi=%f",phi*180/SF_PI);*/
/*    vout("jj",&jj);*/

    nj = rot3d(&jj,nn,tht);
    return nj;
}

/*------------------------------------------------------------*/
vc3d vcp3dfix(vc3d* U, vc3d* V)
/*< vector product of 3D vectors >*/
{
    vc3d W;

    W.dx=(U->dy*V->dz) - (V->dy*U->dz);
    W.dy=(U->dz*V->dx) - (V->dz*U->dx);
    W.dz=(U->dx*V->dy) - (V->dx*U->dy);

    return W;
}

