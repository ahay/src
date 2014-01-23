/*************************************************************************
 * Calculate divergence, curl, polarization operators as well as
 * projection deviation operators in wavenumber domain
 * 1) define wave-vector and devergence operators;
 * 2) solve kelvin-Christoffel equation for the eigenvalue & eigeinvectors;
 * 3) calculate P-waves' polarization operators;
 * 4) calculate P-wave's projection deviation operator.
 *
 *************************************************************************/
/*
   Copyright (C) 2012 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng
     
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
#include "engein2dvti.h"

void divpoldevttip(float **adx,float **adz, float **apx,float **apz, float **apvx,float **apvz,
             float *kx,float *kz, float *kkx,float *kkz, float **taper,          
             int hnkx, int hnkz, float dkx, float dkz,
             double vp2, double vs2, double ep2, double de2, double the, int itaper)
/*< divpoldevttip: divergence, P-wave's polarization operators and projection
  deviation operators for P-wave in VTI media >*/
{
        int   i, j, ii, jj, ik, jk, num;
        double kx2, kz2, k2, rk, sinx, cosx;
        double coss, sins;
        double sum, kxx, kzz;

        double ve[2][2],va[2];

        coss=cos(the);
        sins=sin(the);

        for( i=-hnkx; i<=hnkx ; i++ )
        {
           ik=i+hnkx;

           for( j=-hnkz; j<=hnkz ; j++)
           {
                jk=j+hnkz;

                if(i==0 && j==0 )
                {
                  adx[ik][jk]=0.0;
                  adz[ik][jk]=0.0;
                  apx[ik][jk]=0.0;
                  apz[ik][jk]=0.0;
		  apvx[ik][jk]=1.0;
		  apvz[ik][jk]=1.0;
                  continue;
                }
                // rotatiing according to tilted symmetry axis
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

                if(kzz==0.0) { /*move this sample from zero for cintinuous spectrum */
                  cosx = 0.0001*dkz/rk;
                  sinx = SGN(kxx)*sqrt(1.0-cosx*cosx);
                }

                adx[ik][jk]=sinx;
                adz[ik][jk]=cosx;

                engein2dvti2(ve, va, sinx, cosx, vp2, vs2, ep2, de2);
                
                apx[ik][jk]=(float)ve[0][0];
                apz[ik][jk]=(float)ve[0][1];

		apvx[ik][jk]=(float)(ve[0][0]/sinx);
		apvz[ik][jk]=(float)(ve[0][1]/cosx);
                /*
                if(i!=0)
		  apvx[ik][jk]=(float)(ve[0][0]/sinx);
                else
		  apvx[ik][jk]=1.0;
			
                if(j!=0)
		  apvz[ik][jk]=(float)(ve[0][1]/cosx);
                else
		  apvz[ik][jk]=1.0;
                */
                if(itaper==1){
                  adx[ik][jk] *=taper[ik][jk];
                  adz[ik][jk] *=taper[jk][ik];
                  apx[ik][jk] *=taper[ik][jk];
                  apz[ik][jk] *=taper[jk][ik];

		  apvx[ik][jk] *=taper[ik][jk];
		  apvz[ik][jk] *=taper[jk][ik];
                }
          } /* j loop */
      } /*i loop */

     int nkx, nkz;
     nkx=2*hnkx+1;
     nkz=2*hnkz+1;

     /* interpolating */

     for( i=0; i<nkx; i++ )
     {
          adz[i][hnkz]=(adz[i][hnkz+1] + adz[i][hnkz-1])/2.0;
          adx[i][hnkz]=(adx[i][hnkz+1] + adx[i][hnkz-1])/2.0;
          apz[i][hnkz]=(apz[i][hnkz+1] + apz[i][hnkz-1])/2.0;
          apx[i][hnkz]=(apx[i][hnkz+1] + apx[i][hnkz-1])/2.0;
          apvz[i][hnkz]=(apvz[i][hnkz+1] + apvz[i][hnkz-1])/2.0;
          apvx[i][hnkz]=(apvx[i][hnkz+1] + apvx[i][hnkz-1])/2.0;
     }
     for( j=0; j<nkz ; j++)
     {
          adx[hnkx][j]=(adx[hnkx+1][j] + adx[hnkx-1][j])/2.0;
          adz[hnkx][j]=(adz[hnkx+1][j] + adz[hnkx-1][j])/2.0;
          apx[hnkx][j]=(apx[hnkx+1][j] + apx[hnkx-1][j])/2.0;
          apz[hnkx][j]=(apz[hnkx+1][j] + apz[hnkx-1][j])/2.0;
          apvx[hnkx][j]=(apvx[hnkx+1][j] + apvx[hnkx-1][j])/2.0;
          apvz[hnkx][j]=(apvz[hnkx+1][j] + apvz[hnkx-1][j])/2.0;
     }

       if(fabs(fabs(the)-SF_PI/4.0)<0.001)
       {
         for( i=0; i<nkx; i++)
         for( j=0; j<nkz; j++)
         {
            if(fabs(apvx[i][j])<0.01)
            {
               sum=0.0;
               num=0;
               for(ii=-2;ii<=2;ii++)
               {
                  if(i+ii>=0&&i+ii<nkx)
                     for(jj=-2;jj<=2;jj++)
                     {
                         if(j+jj>=0&&j+jj<nkz)
                         {
                            sum+=apvx[i+ii][j+jj];
                            num++;
                         }
                     }
               }
               apvx[i][j]=sum/num;
            } // end if
            if(fabs(apvz[i][j])<0.01)
            {
               sum=0.0;
               num=0;
               for(ii=-2;ii<=2;ii++)
               {
                  if(i+ii>=0&&i+ii<nkx)
                     for(jj=-2;jj<=2;jj++)
                     {
                         if(j+jj>=0&&j+jj<nkz)
                         {
                            sum+=apvz[i+ii][j+jj];
                            num++;
                         }
                      }
               }
               apvz[i][j]=sum/num;
            } // end if
          }//j loop
       }
}
