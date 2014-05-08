/*************************************************************************
 * Calculate projection deviation operators in wavenumber domain 
 * 1) define wave-vector and devergence operators;
 * 2) solve kelvin-Christoffel equation for the eigenvalue & eigeinvectors
 * 3) calculate projection deviation due to difference between 
 *    wave-vector and P-wave's polarization vector
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
#include "smth2d.h"
#include "engein2dvti.h"

void devvtip(float **apvx,float **apvz, float *kx, float *kz, float *kkx, float*kkz,
             float *kx2, float *kz2, float **taper, int hnkx, int hnkz,
             float dkx, float dkz,
             double vp2, double vs2, double ep2, double de2, int itaper)
/*< devvtip: projection deviation operators for P-wave in VTI media >*/
{
        int   i, j, ik, jk;

        double k2, rk, sinx, cosx;

        double   ve[2][2], va[2];  /*eigeinvector and eigeinvalues*/

	int nkx, nkz;

        for( i=-hnkx; i<=hnkx ; i++ )
        {
           ik=i+hnkx;
           for( j=-hnkz; j<=hnkz ; j++)
           {
                jk=j+hnkz;

                if(i==0&&j==0)
                {
                  apvx[ik][jk]=1.0;
                  apvz[ik][jk]=1.0;
                  continue;
                }
                k2=kx2[ik]+kz2[jk];
                rk=sqrt(k2);

                sinx=kx[ik]/rk;
                cosx=kz[jk]/rk;

                if(kx[ik]==0.0)  { /*move this sample from zero for cintinuous spectrum */
                  sinx = 0.0001*dkx/rk;
                  cosx = SGN(kz[jk])*sqrt(1.0-sinx*sinx);
                }

                if(kz[jk]==0.0) { /*move this sample from zero for cintinuous spectrum */
                  cosx = 0.0001*dkz/rk;    
                  sinx = SGN(kx[ik])*sqrt(1.0-cosx*cosx);
                }

		/*
                //engein2dvti1(ve, va, sinx, cosx, vp2, vs2, ep2, de2, f);
                //sf_warning("Dellinger: va[0]=%f va[1]=%f",va[0],va[1]);
                //sf_warning("ve[0][0]=%f ve[0][1]=%f",ve[0][0],ve[0][1]);
                //sf_warning("ve[1][0]=%f ve[1][1]=%f",ve[1][0],ve[1][1]);
		*/
                
                engein2dvti2(ve, va, sinx, cosx, vp2, vs2, ep2, de2);
		/*
                //sf_warning("2*2Matrix: va[0]=%f va[1]=%f",va[0],va[1]);
                //sf_warning("ve[0][0]=%f ve[0][1]=%f",ve[0][0],ve[0][1]);
                //sf_warning("ve[1][0]=%f ve[1][1]=%f",ve[1][0],ve[1][1]);
                
                //engein2dvti3(ve, va, sinx, cosx, vp2, vs2, ep2, de2);
                //sf_warning("Lapack: va[0]=%f va[1]=%f",va[0],va[1]);
                //sf_warning("ve[0][0]=%f ve[0][1]=%f",ve[0][0],ve[0][1]);
                //sf_warning("ve[1][0]=%f ve[1][1]=%f",ve[1][0],ve[1][1]);
		*/

		apvx[ik][jk]=(float)(ve[0][0]/sinx);
		apvz[ik][jk]=(float)(ve[0][1]/cosx);

                /*
                if(i!=0)
		  apvx[ik][jk]=ve[0][0]/sinx;
                else
		  apvx[ik][jk]=1.0;
			
                if(j!=0)
		  apvz[ik][jk]=ve[0][1]/cosx;
                else
		  apvz[ik][jk]=1.0;
                */
                if(itaper==1){
		  apvx[ik][jk] *=taper[ik][jk];
		  apvz[ik][jk] *=taper[jk][ik];
                }

          } /* j loop */
      } /*i loop */

     nkx=2*hnkx+1;
     nkz=2*hnkz+1;

     /* interpolating */

     for( i=0; i<nkx; i++ )
     {
          apvx[i][hnkz]=(apvx[i][hnkz+1] + apvx[i][hnkz-1])/2.0;
          apvz[i][hnkz]=(apvz[i][hnkz+1] + apvz[i][hnkz-1])/2.0;
     }
     for( j=0; j<nkz ; j++)
     {
          apvx[hnkx][j]=(apvx[hnkx+1][j] + apvx[hnkx-1][j])/2.0;
          apvz[hnkx][j]=(apvz[hnkx+1][j] + apvz[hnkx-1][j])/2.0;
     }
}

void devvtips(float **apvx,float **apvz, float *kx, float *kz, float *kkx, float*kkz,
             float *kx2, float *kz2, float **taper, int hnkx, int hnkz,
             float dkx, float dkz,
             double vp2, double vs2, double ep2, double de2, int itaper)
/*< devvtips: projection deviation operators for SV-wave in VTI media >*/
{
        int   i, j, ik, jk;

        double k2, rk, sinx, cosx;

        double   ve[2][2], va[2];  /*eigeinvector and eigeinvalues*/

	int nkx, nkz;

        for( i=-hnkx; i<=hnkx ; i++ )
        {
           ik=i+hnkx;
           for( j=-hnkz; j<=hnkz ; j++)
           {
                jk=j+hnkz;

                if(i==0&&j==0)
                {
                  apvx[ik][jk]=1.0;
                  apvz[ik][jk]=1.0;
                  continue;
                }
                k2=kx2[ik]+kz2[jk];
                rk=sqrt(k2);

                sinx=kx[ik]/rk;
                cosx=kz[jk]/rk;

                if(kx[ik]==0.0)  { /*move this sample from zero for cintinuous spectrum */
                  sinx = 0.01*dkx/rk;
                  cosx = SGN(kz[jk])*sqrt(1.0-sinx*sinx);
                }

                if(kz[jk]==0.0) { /*move this sample from zero for cintinuous spectrum */
                  cosx = 0.01*dkz/rk;    
                  sinx = SGN(kx[ik])*sqrt(1.0-cosx*cosx);
                }

                engein2dvti2(ve, va, sinx, cosx, vp2, vs2, ep2, de2);

		apvx[ik][jk]=(float)(ve[1][0]/sinx);
		apvz[ik][jk]=(float)(ve[1][1]/cosx);

                /*
                if(i!=0)
		  apvx[ik][jk]=ve[0][0]/sinx;
                else
		  apvx[ik][jk]=1.0;
			
                if(j!=0)
		  apvz[ik][jk]=ve[0][1]/cosx;
                else
		  apvz[ik][jk]=1.0;
                */
                if(itaper==1){
		  apvx[ik][jk] *=taper[ik][jk];
		  apvz[ik][jk] *=taper[jk][ik];
                }

          } /* j loop */
      } /*i loop */

     nkx=2*hnkx+1;
     nkz=2*hnkz+1;

     /* interpolating */

     for( i=0; i<nkx; i++ )
     {
          apvx[i][hnkz]=(apvx[i][hnkz+1] + apvx[i][hnkz-1])/2.0;
          apvz[i][hnkz]=(apvz[i][hnkz+1] + apvz[i][hnkz-1])/2.0;
     }
     for( j=0; j<nkz ; j++)
     {
          apvx[hnkx][j]=(apvx[hnkx+1][j] + apvx[hnkx-1][j])/2.0;
          apvz[hnkx][j]=(apvz[hnkx+1][j] + apvz[hnkx-1][j])/2.0;
     }

}
