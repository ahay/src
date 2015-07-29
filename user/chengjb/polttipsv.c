/*************************************************************************
 * Calculate polarization operators in wavenumber domain 
 * 1) solve kelvin-Christoffel equation for the eigenvalue & eigeinvectors
 * 2) calculate P-wave's polarization operators 
 * 3) calculate SV-wave's polarization operators 
 *
 *************************************************************************/
/*
   Copyright (C) 2012 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng
   2012.3.6
   Modified on Madagascar at 2012.8.1
     
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
#include "smth2d.h"
#include "engein2dvti.h"

void polttipsv(float **apx,float **apz, float **apxs,float **apzs,
             float *kx,float *kz, float *kkx,float *kkz,
             float **taper, int hnkx, int hnkz, float dkx, float dkz,
             double vp2, double vs2, double ep2, double de2, double the, int itaper)
/*< polvtipsv: P- and SV-wave polarization operators in VTI media >*/
{
        int   i, j, ik, jk;
        double k2, rk, sinx, cosx, coss, sins;
        double kxx, kzz, kx2, kz2;
        double ve[2][2], va[2];
	int nkx, nkz;

        coss=cos(the);
        sins=sin(the);

        for( i=-hnkx; i<=hnkx ; i++ )
        {
           ik=i+hnkx;
           for( j=-hnkz; j<=hnkz ; j++)
           {
                jk=j+hnkz;
                if(i==0&&j==0)
                {
                  apx[ik][jk]=0.0;
                  apz[ik][jk]=0.0;
                  apxs[ik][jk]=0.0;
                  apzs[ik][jk]=0.0;
                  continue;
                }

                /* rotatiing according to tilted symmetry axis */
                kxx=kx[ik]*coss+kz[jk]*sins;
                kzz=kz[jk]*coss-kx[ik]*sins;
                kx2=kxx*kxx;
                kz2=kzz*kzz;
                k2=kx2+kz2;
                rk=sqrt(k2);

                sinx=kxx/rk;
                cosx=kzz/rk;

                if(kxx==0.0)  { /*move this sample from zero for cintinuous spectrum */
                  sinx = 0.0001*dkx/rk;
                  cosx = SGN(kzz)*sqrt(1.0-sinx*sinx);
                }

                if(kz[jk]==0.0) { /*move this sample from zero for cintinuous spectrum */
                  cosx = 0.0001*dkz/rk;
                  sinx = SGN(kxx)*sqrt(1.0-cosx*cosx);
                }

                engein2dvti2(ve, va, sinx, cosx, vp2, vs2, ep2, de2);

                apx[ik][jk]=(float)ve[0][0];
                apz[ik][jk]=(float)ve[0][1];
                if(ve[0][0]*sinx + ve[0][1]*cosx <=0.0)
                {
                     apx[ik][jk]= (float)(-ve[0][0]);
                     apz[ik][jk]= (float)(-ve[0][1]);
                }
                apx[ik][jk] *= taper[ik][jk];
                apz[ik][jk] *= taper[jk][ik];

                apxs[ik][jk]=(float)ve[1][0];
                apzs[ik][jk]=(float)ve[1][1];
                if(ve[1][0]*cosx - ve[1][1]*sinx <=0.0)
                {
                     apxs[ik][jk]= (float)(-ve[1][0]);
                     apzs[ik][jk]= (float)(-ve[1][1]);
                }
                if(itaper==1){
                  apxs[ik][jk] *= taper[jk][ik];
                  apzs[ik][jk] *= taper[ik][jk];
                }
          } /* j loop */
      } /*i loop */

     nkx=2*hnkx+1;
     nkz=2*hnkz+1;
     /* interpolating */
     for( i=0; i<nkx; i++ )
     {
          apx[i][hnkz]=(apx[i][hnkz+1] + apx[i][hnkz-1])/2.0;
          apz[i][hnkz]=(apz[i][hnkz+1] + apz[i][hnkz-1])/2.0;
          apxs[i][hnkz]=(apxs[i][hnkz+1] + apxs[i][hnkz-1])/2.0;
          apzs[i][hnkz]=(apzs[i][hnkz+1] + apzs[i][hnkz-1])/2.0;
     }
     for( j=0; j<nkz ; j++)
     {
          apx[hnkx][j]=(apx[hnkx+1][j] + apx[hnkx-1][j])/2.0;
          apz[hnkx][j]=(apz[hnkx+1][j] + apz[hnkx-1][j])/2.0;
          apxs[hnkx][j]=(apxs[hnkx+1][j] + apxs[hnkx-1][j])/2.0;
          apzs[hnkx][j]=(apzs[hnkx+1][j] + apzs[hnkx-1][j])/2.0;
     }
}
