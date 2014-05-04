/*************************************************************************
 * Calculate eigeinvalues and eigeinvectors for 2D VTI media
 * Method 1: analytic solution based on Dellinger's expression
 * Method 2: using 2*2 matrix analytic solution
 * Method 3: using general N*N matrix Lapack solution
 *
 * Note:  Method 2 & 3 are similar, they solve same equation using different codes
 *
 *************************************************************************/
/*
   Copyright (C) 2012 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng, Wei Kang and Tengfei Wang
     
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
#include "eigen2x2.h"

void engein2dvti1(double ve[2][2], double va[2], double sinx, double cosx, double vp2, double vs2,
                 double ep2, double de2, double f)
/*< engein2dvti1: Calculate eigeinvalues and eigeinvectors for 2D VTI media
                  using analytic solution based on Dellinger's expression >*/
{
        double sin2, cos2;
        double d33, d11, psi2, psi, d33d11, d33d11_2, sin2cos2, sin2cos2_4, tmpa,tmpb,u1,u2,u1u2;

        sin2=sinx*sinx;
        cos2=1.0-sin2;
/*        cos2a=1.0-2*sin2; */
/*        sin4=sin2*sin2;   */
        sin2cos2=sin2*cos2;
        sin2cos2_4=4*sin2cos2;

        /* Dellinger's direct method PhD Chapter2 P.12 */ 
        // d33 = C33 - C55;
        // d11 = C11 - C55;
        // psi = C13 + C55
        d33=vp2-vs2;
        d11=ep2*vp2-vs2;
        psi2=d33*(de2*vp2-vs2);
        psi=sqrtf(psi2);
        d33d11=d33*cos2-d11*sin2;
        d33d11_2=d33d11*d33d11;
        tmpa=vs2+vp2*cos2+(ep2*vp2)*sin2;
        tmpb=sqrtf(d33d11_2+sin2cos2_4*psi2);

        va[0]=sqrtf(0.5*(tmpa+tmpb));  // P-wave phase vlocity
        va[1]=sqrtf(0.5*(tmpa-tmpb));  // SV-wave phase vlocity

        u1=2*psi*sqrtf(sin2*cos2);
        u2=sqrtf(d33d11_2+sin2cos2_4*psi2) + d33d11;

        /* normalize the polarization vector */
        u1u2=sqrt(u1*u1+u2*u2);
        if(u1u2==0)
        {
          u1=0.0;
          u2=0.0;
        }else
        {
          u1=u1/u1u2;
          u2=u2/u1u2;
        }
        /* get the closest direction to k */
        if(u1*sinx + u2*cosx <0)
        {
           u2 = -u2;
           u1 = -u1;
        }
        ve[0][0]=u1;
        ve[0][1]=u2;

        u1=ve[0][1];
        u2=-ve[0][0];

        /* get the closest direction to k */
        if(u1*cosx - u2*sinx <0)
        {
           u2 = -u2;
           u1 = -u1;
        }
        ve[1][0]=u1;
        ve[1][1]=u2;
}

void engein2dvti2(double ve[2][2], double va[2], double kx, double kz, double vp2, double vs2,
                 double ep2, double de2)
/*< engein2dvti2: Calculate eigeinvalues and eigeinvectors for 2D VTI media
                  using 2*2 matrix analytic solution>*/
{
        double u1, u2, c11, c33, c44, c13c44, a11, a12, a22;
        double a[2][2];

        c33=vp2;
        c44=vs2;
        c11=ep2*c33;
        //c13=sqrt(2*c33*(c33-c44)*de+(c33-c44)*(c33-c44))-c44;
        //c13=sqrt((2*c33*de+(c33-c44))*(c33-c44))-c44;
        //c13=sqrt((2*c33*de+c33-c44)*(c33-c44))-c44;
        //c13=sqrt((2*de+1.0)*c33-c44)*(c33-c44))-c44;
        c13c44=sqrt((de2*c33-c44)*(c33-c44));

        a11= c11*kx*kx+c44*kz*kz;
        a12= c13c44*kx*kz;
        a22= c44*kx*kx+c33*kz*kz;

        a[0][0] = a11;
        a[0][1] = a12;
        a[1][0] = a12;
        a[1][1] = a22;

        dsolveSymmetric22(a, ve, va);

        u1=ve[0][0];
        u2=ve[0][1];

        /* get the closest direction to k */
        if(u1*kx + u2*kz <0) {
           u2 = -u2;
           u1 = -u1;
        }
        ve[0][0]=u1;
        ve[0][1]=u2;

        u1=ve[1][0];
        u2=ve[1][1];

        /* get the closest direction to k */
        if(u1*kz - u2*kx <0) {
           u2 = -u2;
           u1 = -u1;
        }
        ve[1][0]=u1;
        ve[1][1]=u2;
}
