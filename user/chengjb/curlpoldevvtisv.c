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

void curlpoldevvtisv(float **acx,float **acz, float **asx,float **asz,float **asvx,float **asvz,
                   float *kx,float *kz, float *kkx,float *kkz,
                   float *kx2, float *kz2, float **taper, int hnkx, int hnkz, float dkx, float dkz,
                   double vp2, double vs2, double ep2, double de2, int itaper)
/*< curlpoldevvtip: curl, qSV-wave's polarization operators and projection
  deviation operators for qSV-wave in VTI media >*/
{
        int   i, j, ik, jk;

        double k2, rk, sinx, cosx;
        double ve[2][2], va[2];  /*eigeinvector and eigeinvalues*/

        for( i=-hnkx; i<=hnkx ; i++ )
        {
           ik=i+hnkx;
           for( j=-hnkz; j<=hnkz ; j++)
           {
                jk=j+hnkz;

                if(i==0&&j==0)
                {
                  acx[ik][jk]=0.0;
                  acz[ik][jk]=0.0;
                  asx[ik][jk]=0.0;
                  asz[ik][jk]=0.0;
                  asvx[ik][jk]=1.0;
                  asvz[ik][jk]=1.0;
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

                acx[ik][jk]=cosx;
                acz[ik][jk]=-sinx;
                //sf_warning("acx=%f acz=%f",acx[ik][jk],acz[ik][jk]);

                //engein2dvti1(ve, va, sinx, cosx, vp2, vs2, ep2, de2, f);
                //sf_warning("Dellinger: va[0]=%f va[1]=%f",va[0],va[1]);
                //sf_warning("ve[0][0]=%f ve[0][1]=%f",ve[0][0],ve[0][1]);
                //sf_warning("ve[1][0]=%f ve[1][1]=%f",ve[1][0],ve[1][1]);
                
                engein2dvti2(ve, va, sinx, cosx, vp2, vs2, ep2, de2);
                //sf_warning("2*2Matrix: va[0]=%f va[1]=%f",va[0],va[1]);
                //sf_warning("ve[0][0]=%f ve[0][1]=%f",ve[0][0],ve[0][1]);
                //sf_warning("ve[1][0]=%f ve[1][1]=%f",ve[1][0],ve[1][1]);
                
                //engein2dvti3(ve, va, sinx, cosx, vp2, vs2, ep2, de2);
                //sf_warning("Lapack: va[0]=%f va[1]=%f",va[0],va[1]);
                //sf_warning("ve[0][0]=%f ve[0][1]=%f",ve[0][0],ve[0][1]);
                //sf_warning("ve[1][0]=%f ve[1][1]=%f",ve[1][0],ve[1][1]);

                asx[ik][jk]=(float)ve[1][0];
                asz[ik][jk]=(float)ve[1][1];
                //sf_warning("apx=%f apz=%f",(float)ve[0][0],(float)ve[0][1]);
                //sf_warning("asx=%f asz=%f",asx[ik][jk],asz[ik][jk]);

				asvx[ik][jk]=(float)(ve[1][0]/cosx);
				asvz[ik][jk]=(float)(-1*ve[1][1]/sinx);
				//asvx[ik][jk]=(float)(-1*ve[1][0]/cosx);
				//asvz[ik][jk]=(float)(ve[1][1]/sinx);
                //sf_warning("asvx=%f asvz=%f",asvx[ik][jk],asvz[ik][jk]);
                /*
                if(i!=0)
		  asvx[ik][jk]=(float)(ve[0][0]/sinx);
                else
		  asvx[ik][jk]=1.0;
			
                if(j!=0)
		  asvz[ik][jk]=(float)(ve[0][1]/cosx);
                else
		  asvz[ik][jk]=1.0;
                */
                if(itaper==1){
				  acx[ik][jk] *=taper[ik][jk];
				  acz[ik][jk] *=taper[jk][ik];
				  asx[ik][jk] *=taper[ik][jk];
				  asz[ik][jk] *=taper[jk][ik];
				  asvx[ik][jk] *=taper[ik][jk];
				  asvz[ik][jk] *=taper[jk][ik];
                }

          } /* j loop */
      } /*i loop */

     /* interpolating */

     int nkx, nkz;
     nkx=2*hnkx+1;
     nkz=2*hnkz+1;
     for( i=0; i<nkx; i++ )
     {
          acx[i][hnkz]=(acx[i][hnkz+1] + acx[i][hnkz-1])/2.0;
          acz[i][hnkz]=(acz[i][hnkz+1] + acz[i][hnkz-1])/2.0;
          asx[i][hnkz]=(asx[i][hnkz+1] + asx[i][hnkz-1])/2.0;
          asz[i][hnkz]=(asz[i][hnkz+1] + asz[i][hnkz-1])/2.0;
          asvx[i][hnkz]=(asvx[i][hnkz+1] + asvx[i][hnkz-1])/2.0;
          asvz[i][hnkz]=(asvz[i][hnkz+1] + asvz[i][hnkz-1])/2.0;
     }
     for( j=0; j<nkz ; j++)
     {
          acx[hnkx][j]=(acx[hnkx+1][j] + acx[hnkx-1][j])/2.0;
          acz[hnkx][j]=(acz[hnkx+1][j] + acz[hnkx-1][j])/2.0;
          asx[hnkx][j]=(asx[hnkx+1][j] + asx[hnkx-1][j])/2.0;
          asz[hnkx][j]=(asz[hnkx+1][j] + asz[hnkx-1][j])/2.0;
          asvx[hnkx][j]=(asvx[hnkx+1][j] + asvx[hnkx-1][j])/2.0;
          asvz[hnkx][j]=(asvz[hnkx+1][j] + asvz[hnkx-1][j])/2.0;
     }

}
