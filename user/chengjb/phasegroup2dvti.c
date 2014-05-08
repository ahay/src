/*************************************************************************
 * Calculate phase velocity, group angle amd group velocity for
 * qP, qSV-waves in 2D VTI media
 *************************************************************************/
/*
 
  Copyright (C) 2012 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng, Wei Kang and Tengfei Wang.
     
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

#include "_cjb.h"

float vpphase2dvti(float vp0, float vs0, float eps, float del, float ap)
/*< vpphase2dvti: calculate phase velocity for given phase angle of qP-wave in 2D VTI media >*/
{
        float alpha = vs0/vp0;
        float f = 1-alpha*alpha;
        float vpp;

        vpp = vp0*sqrt(1+eps*sin(ap)*sin(ap)-f/2.0
           + f/2.0*sqrt(pow(1+2*eps*sin(ap)*sin(ap)/f,2)-2*(eps-del)*sin(2*ap)*sin(2*ap)/f));

	return vpp;
}

float vsphase2dvti(float vp0, float vs0, float eps, float del, float ap)
/*< vsphase2dvti: calculate phase velocity for given phase angle of qSV-wave in 2D VTI media >*/
{
        float alpha = vs0/vp0;
        float f = 1-alpha*alpha;
        float vsp;

        vsp = vp0*sqrt(1+eps*sin(ap)*sin(ap)-f/2.0
           - f/2.0*sqrt(pow(1+2*eps*sin(ap)*sin(ap)/f,2)-2*(eps-del)*sin(2*ap)*sin(2*ap)/f));

	return vsp;
}

float vpgroup2dvti(float vp0, float vs0, float eps, float del, float ap)
/*< vpgroup2dvti: calculate group velocity from phase velocity of qP-wave in 2D VTI media >*/
{
        float alpha = vs0/vp0;
        float f = 1-alpha*alpha;

	float vpp = vpphase2dvti(vp0, vs0, eps, del, ap);
        
	float A = eps*sin(2*ap)
                + f*f/4.0/sqrt(4*eps*eps*pow(sin(ap),4)+f*4.0*sin(ap)*sin(ap)*(2*del*cos(ap)*cos(ap)-eps*cos(2*ap))+f*f)
                * (4*eps*sin(2*ap)/f+8*eps*eps*sin(ap)*sin(ap)*sin(2*ap)/(f*f)-8*(eps-del)*sin(2*ap)*cos(2*ap)/f);

        float dvda = A*vp0*vp0/vpp/2.0;

        float vpg = vpp*sqrt(1+pow(dvda/vpp,2));
        return vpg;
}

float vsgroup2dvti(float vp0, float vs0, float eps, float del, float ap)
/*< vsgroup2dvti: calculate group velocity from phase velocity of qSV-wave in 2D VTI media >*/
{
        float alpha = vs0/vp0;
        float f = 1-alpha*alpha;

	float vsp = vsphase2dvti(vp0, vs0, eps, del, ap);
        
        float A = eps*sin(2*ap)
                - f*f/4.0/sqrt(4*eps*eps*pow(sin(ap),4)+f*4.0*sin(ap)*sin(ap)*(2*del*cos(ap)*cos(ap)-eps*cos(2*ap))+f*f)
                * (4*eps*sin(2*ap)/f+8*eps*eps*sin(ap)*sin(ap)*sin(2*ap)/(f*f)-8*(eps-del)*sin(2*ap)*cos(2*ap)/f);

        float dvda = A*vp0*vp0/vsp/2.0;

        float vsg=vsp*sqrt(1+pow(dvda/vsp,2));
        return vsg;
}

float apgroup2dvti(float vp0, float vs0, float eps, float del, float ap)
/*< apgroup2dvti: calculate group angle from phase velocity and angle of qP-wave in 2D VTI media >*/
{
        float alpha = vs0/vp0;
        float f = 1-alpha*alpha;

	float vpp = vpphase2dvti(vp0, vs0, eps, del, ap);

	float A = eps*sin(2*ap)
                + f*f/4.0/sqrt(4*eps*eps*pow(sin(ap),4)+f*4.0*sin(ap)*sin(ap)*(2*del*cos(ap)*cos(ap)-eps*cos(2*ap))+f*f)
                * (4*eps*sin(2*ap)/f+8*eps*eps*sin(ap)*sin(ap)*sin(2*ap)/(f*f)-8*(eps-del)*sin(2*ap)*cos(2*ap)/f);

        float dvda = A*vp0*vp0/vpp/2.0;

	float B = 1-tan(ap)*dvda/vpp;
	float D = tan(ap)+dvda/vpp;

	float C, apg;

	if(fabs(B)<1e-30)
		B=1e-30;

	C=D/B;
	if(fabs(C)>1e30)
		C=1e30;

	apg=atan(C);
	if(apg<0)
		apg+=SF_PI;
	
	return apg;
}

float asgroup2dvti(float vp0, float vs0, float eps, float del, float ap)
/*< asgroup2dvti: calculate group angle from phase velocity and angle of qSV-wave in 2D VTI media >*/
{
        float alpha = vs0/vp0;
        float f = 1-alpha*alpha;

	float vsp = vsphase2dvti(vp0, vs0, eps, del, ap);

        float A = eps*sin(2*ap)
                - f*f/4.0/sqrt(4*eps*eps*pow(sin(ap),4)+f*4.0*sin(ap)*sin(ap)*(2*del*cos(ap)*cos(ap)-eps*cos(2*ap))+f*f)
                * (4*eps*sin(2*ap)/f+8*eps*eps*sin(ap)*sin(ap)*sin(2*ap)/(f*f)-8*(eps-del)*sin(2*ap)*cos(2*ap)/f);

        float dvda = A*vp0*vp0/vsp/2.0;

	float B=1-tan(ap)*dvda/vsp;
	float D=tan(ap)+dvda/vsp;

	float C, asg;

	if(fabs(B)<1e-20)
		B=1e-20;

	C=D/B;
	if(fabs(C)>1e5)
		C=1e5;

	asg=atan(C);
	if(asg<0)
		asg+=SF_PI;
	
	return asg;
}

void vapgroup2dvti(float vp0, float vs0, float eps, float del, float ap, float *vpg, float *apg)
/*< vapgroup2dvti: calculate group velocity and angle from phase velocity for
 * given phase angle of qP-wave in 2D VTI media >*/
{
        float alpha = vs0/vp0;
        float f = 1-alpha*alpha;

	float vpp = vpphase2dvti(vp0, vs0, eps, del, ap);
        
	float A = eps*sin(2*ap)
                + f*f/4.0/sqrt(4*eps*eps*pow(sin(ap),4)+f*4.0*sin(ap)*sin(ap)*(2*del*cos(ap)*cos(ap)-eps*cos(2*ap))+f*f)
                * (4*eps*sin(2*ap)/f+8*eps*eps*sin(ap)*sin(ap)*sin(2*ap)/(f*f)-8*(eps-del)*sin(2*ap)*cos(2*ap)/f);

        float dvda = A*vp0*vp0/vpp/2.0;

	float B, D, C;

        *vpg = vpp*sqrt(1+pow(dvda/vpp,2));

	B=1-tan(ap)*dvda/vpp;
	D=tan(ap)+dvda/vpp;

	if(fabs(B)<1e-30)
		B=1e-30;

	C=D/B;

	if(fabs(C)>1e30)
		C=1e30;

	*apg=atan(C);
	if(*apg<0)
		*apg+=SF_PI;
}

void vasgroup2dvti(float vp0, float vs0, float eps, float del, float ap, float *vsg, float *asg)
/*< vasgroup2dvti: calculate group velocity and angle from phase velocity for
 * given phase angle of qSV-wave in 2D VTI media >*/
{
        float alpha = vs0/vp0;
        float f = 1-alpha*alpha;

	float vsp = vsphase2dvti(vp0, vs0, eps, del, ap);
        
        float A = eps*sin(2*ap)
                - f*f/4.0/sqrt(4*eps*eps*pow(sin(ap),4)+f*4.0*sin(ap)*sin(ap)*(2*del*cos(ap)*cos(ap)-eps*cos(2*ap))+f*f)
                * (4*eps*sin(2*ap)/f+8*eps*eps*sin(ap)*sin(ap)*sin(2*ap)/(f*f)-8*(eps-del)*sin(2*ap)*cos(2*ap)/f);

        float dvda = A*vp0*vp0/vsp/2.0;

	float B, D, C;

        *vsg = vsp*sqrt(1+pow(dvda/vsp,2));

	B=1-tan(ap)*dvda/vsp;
	D=tan(ap)+dvda/vsp;

	if(fabs(B)<1e-30)
		B=1e-30;

	C=D/B;
	if(fabs(C)>1e20)
		C=1e20;

	*asg=atan(C);
	if(*asg<0)
		*asg+=SF_PI;
}
