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

#include "green.h"

#define EPS 1.0e-6

static float v0,vgrad;

void GreenInit(float v0_in, float vgrad_in)
/*< initialize >*/
{
    v0=v0_in;
    vgrad=vgrad_in;
}

void ZeroArray(float *xx,int n)
/*< zero an array >*/
{
    int i;
    for (i=0; i<n; i++) xx[i]=0.;
    return;
}

float Vel(float x1,float x2,float x3)
{
    float vel;
    vel=v0+vgrad*x3;
    return (vel);
}

sf_complex Green(float r1,float r2,float r3,float s1,float s2,float s3,float omega) 
/*< Green's function >*/
{
    double tt,amp;
    sf_complex val;

    GreenTtAmp(r1,r2,r3,s1,s2,s3,&tt,&amp);
#ifdef SF_HAS_COMPLEX_H
    val=amp*cexpf(sf_cmplx(0.,omega*tt)); 
#else
    val=sf_crmul(cexpf(sf_cmplx(0.,omega*tt)),amp); 
#endif
    return (val);
}

void GreenTtAmp(float r1,float r2,float r3,float s1,float s2,float s3,double *tt,double *amp)
/*< traveltime-amplitude >*/
{
    double x2,z1,z2,beta,xr;
    double thetas,thetar;

    if (vgrad < EPS) {
	xr=sqrt((s1-r1)*(s1-r1) + (s2-r2)*(s2-r2) + (s3-r3)*(s3-r3));
	*tt=xr/v0;
	*amp=1./SF_MAX(xr,EPS);
    } else {

	x2=sqrt((s1-r1)*(s1-r1) + (s2-r2)*(s2-r2));
	if (x2<EPS) { 
	    *tt=fabs(r3-s3)/(0.5*s3*vgrad+0.5*r3*vgrad+v0);
	    *amp=1./((*tt)*v0);
	    return;
	}

	z1=SF_MIN(s3,r3); z1+=v0/vgrad;
	z2=SF_MAX(s3,r3); z2+=v0/vgrad;

	beta=0.5*(z1+z2)/x2;
	xr=0.5*x2-beta*(z1-z2);

	thetas=atan(z1/xr);
	if (fabs(xr-x2)<EPS) thetar=SF_PI/2.;
	else                    thetar=atan(z2/(xr-x2));

	if (xr > x2)
	    *tt= (log(fabs(tan(thetar/2)))-log(fabs(tan(thetas/2))))/vgrad;
	else if (xr < x2)
	    *tt=-(log(fabs(tan(thetar/2)))+log(fabs(tan(thetas/2))))/vgrad;
	else 
	    *tt=-(log(fabs(tan(thetas/2))));

	*amp=1./(v0*SF_MAX(*tt,EPS));
    }
    return;
}

