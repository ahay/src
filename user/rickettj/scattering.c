/*
  Copyright (C) 2000 The Board of Trustees of Stanford University
  
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

#include "scattering.h"
#include "green.h"

static int nt,nw,border;
static float *u,*du,*ut,*utt;
static float ot,dt,ow,dw,sx,sy,sz,rx,ry,rz,px,py,pz,dv;

void ScatteringInit(float dv1,
		    int nt1,float ot1,float dt1,
		    int nw1,float ow1,float dw1,
		    float sx1,float sy1,float sz1,
		    float rx1,float ry1,float rz1,
		    float px1,float py1,float pz1)
/*< initialize >*/
{
    dv=dv1;

    nt=nt1; ot=ot1, dt=dt1;    
    nw=nw1; ow=ow1; dw=dw1;
    sx=sx1; sy=sy1; sz=sz1;
    rx=rx1; ry=ry1; rz=rz1;
    px=px1; py=py1; pz=pz1;

    u  =sf_floatalloc(nt);
    du =sf_floatalloc(nt);
    ut =sf_floatalloc(nt);
    utt=sf_floatalloc(nt);

    border=nt/5;
    sf_warning("Border width: %d",border);
    return;
}

void BornScatteredField(float xx,float xy,float xz,float *output)
/*< Born scattering field >*/
{
    int iw;
    double omega,scale;

    sf_complex  U,dU;
    sf_complex val,fkern;

    ZeroArray(u,  nt); 
    ZeroArray(du, nt); 
    *output=0.;

    for (iw=0; iw<nw; iw++) {
	omega=ow+dw*iw;

	scale =cos((SF_PI/2)*(((double) iw+1)/((double) nw+1)));
	scale*=scale*omega;

#ifdef SF_HAS_COMPLEX_H
	/* background field */
	U=Green(xx,xy,xz,sx,sy,sz,omega)*scale;
	    
	/* scattered field */
	val=Green(px,py,pz,sx,sy,sz,omega)*scale;
	val *= Green(xx,xy,xz,px,py,pz,omega);
	dU = val*(-omega*omega*dv);

	U += dU;
	fkern=cexpf(sf_cmplx(0.,-omega*0.6));
	val=fkern*U;
#else
	/* background field */
	U=sf_crmul(Green(xx,xy,xz,sx,sy,sz,omega),scale);
	    
	/* scattered field */
	val=sf_crmul(Green(px,py,pz,sx,sy,sz,omega),scale);
	val=sf_cmul(val,Green(xx,xy,xz,px,py,pz,omega));
	dU=sf_crmul(val,-omega*omega*dv);

	U=sf_cadd(U,dU);
	fkern=cexpf(sf_cmplx(0.,-omega*0.6));
	val=sf_cmul(fkern,U);
#endif

	*output+=crealf(val);
    }
    return;
}    


void RytovSensitivity(float xx,float xy,float xz,float *output)
/*< Rytov sensitivity >*/
{
    int iw;
    double omega,scale;
    sf_complex U,dU;
    sf_complex val;

    *output=0.;

    for (iw=0; iw<nw; iw++) {
	omega=ow+dw*iw;

	scale =cos((SF_PI/2)*(((double) iw+1)/((double) nw+1)));
	scale*=scale;

#ifdef SF_HAS_COMPLEX_H
	/* background field */
	U=Green(rx,ry,rz,sx,sy,sz,omega)*scale;
	    
	/* scattered field */
	val=Green(xx,xy,xz,sx,sy,sz,omega)*scale;
	val*=Green(rx,ry,rz,xx,xy,xz,omega);
	dU=val*(-omega*omega*dv);

	val=U*omega;
	*output -= scale*cimagf(dU/U);
#else
	/* background field */
	U=sf_crmul(Green(rx,ry,rz,sx,sy,sz,omega),scale);
	    
	/* scattered field */
	val=sf_crmul(Green(xx,xy,xz,sx,sy,sz,omega),scale);
	val=sf_cmul(val,Green(rx,ry,rz,xx,xy,xz,omega));
	dU=sf_crmul(val,-omega*omega*dv);

	val=sf_crmul(U,omega);
	*output -= scale*cimagf(sf_cdiv(dU,U));
#endif
    }
  
    return;
}

void BornSensitivity(float xx,float xy,float xz,float *output)
/*< Born sensitivity >*/
{
    int iw,it;
    float tt,top,bottom;
    double omega,scale;
    sf_complex val,fkern;
    sf_complex  U,dU,Ut,Utt;

    ZeroArray(u,  nt); 
    ZeroArray(du, nt); 
    ZeroArray(ut, nt); 
    ZeroArray(utt,nt); 

    for (iw=0; iw<nw; iw++) {
	omega=ow+dw*iw;

	scale =cos((SF_PI/2)*(((double) iw+1)/((double) nw+1)));
	scale*=scale;

#ifdef SF_HAS_COMPLEX_H
	/* background field */
	U=Green(rx,ry,rz,sx,sy,sz,omega)*scale;
	    
	/* scattered field */
	val=Green(xx,xy,xz,sx,sy,sz,omega)*scale;
	val*=Green(rx,ry,rz,xx,xy,xz,omega);
	dU=val*(-omega*omega*dv);

	/* field derivatives */
	Ut =U*sf_cmplx(0.,omega);
	Utt=U*(-omega*omega);
#else
	/* background field */
	U=sf_crmul(Green(rx,ry,rz,sx,sy,sz,omega),scale);
	    
	/* scattered field */
	val=sf_crmul(Green(xx,xy,xz,sx,sy,sz,omega),scale);
	val=sf_cmul(val,Green(rx,ry,rz,xx,xy,xz,omega));
	dU=sf_crmul(val,-omega*omega*dv);

	/* field derivatives */
	Ut =sf_cmul (U,sf_cmplx(0.,omega));
	Utt=sf_crmul(U,-omega*omega);
#endif

	/* slow ifft */
	for (it=0; it<nt; it++) {
	    tt=it*dt+ot;
	    
	    fkern=cexpf(sf_cmplx(0.,-omega*tt));

#ifdef SF_HAS_COMPLEX_H    
	    val=fkern*U;     u[it]+=crealf(val);
	    val=fkern*dU;   du[it]+=crealf(val);
	    val=fkern*Ut;   ut[it]+=crealf(val);
	    val=fkern*Utt; utt[it]+=crealf(val);
#else
	    val=sf_cmul(fkern,  U);   u[it]+=crealf(val);
	    val=sf_cmul(fkern, dU);  du[it]+=crealf(val);
	    val=sf_cmul(fkern, Ut);  ut[it]+=crealf(val);
	    val=sf_cmul(fkern,Utt); utt[it]+=crealf(val);
#endif
	}
    }

    top=0; bottom=0.;
    for (it=0; it<nt; it++) {
	top    += ut [it]*du[it]*dt;
	bottom += utt[it]* u[it]*dt;
    }
    *output=top/bottom;
  
    return;
}





