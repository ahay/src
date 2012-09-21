/*************************************************************************
 * Calculate divergence, curl, polarization operators as well as
 * projection deviation operators in wavenumber domain
 * 1) define wave-vector and devergence operators;
 * 2) solve kelvin-Christoffel equation for the eigenvalue & eigeinvectors;
 * 3) calculate P-waves' polarization operators;
 * 4) calculate P-wave's projection deviation operator.
 *
 *    Copyright: Tongji University (Jiubing Cheng)
 *    2012.3.2
 *************************************************************************/
#include <rsf.h>

#include "_cjb.h"
#include "smth2d.h"
#include "engein2dvti.h"

void divpoldevvtip(float **adx,float **adz, float **apx,float **apz,float **apvx,float **apvz,
                   float *kx,float *kz, float *kkx,float *kkz,
                   float *kx2, float *kz2, float **taper, int hnkx, int hnkz,
                   float vp2, float vs2, float ep2, float de2, int itaper)
/*< divpoldevvtip: divergence, P-wave's polarization operators and projection
  deviation operators for P-wave in VTI media >*/
{
        int   i, j, ik, jk;
        float k2, rk, sinx, cosx;

        float   ve[2][2], va[2];  /*eigeinvector and eigeinvalues*/

        //float f=1.0-vs2/vp2;

        for( i=-hnkx; i<=hnkx ; i++ )
        {
           ik=i+hnkx;
           for( j=-hnkz; j<=hnkz ; j++)
           {
                jk=j+hnkz;

                if(i==0&&j==0)
                {
                  adx[ik][jk]=0.0;
                  adz[ik][jk]=0.0;
                  apx[ik][jk]=0.0;
                  apz[ik][jk]=0.0;
                  apvx[ik][jk]=1.0;
                  apvz[ik][jk]=1.0;
                  continue;
                }
                k2=kx2[ik]+kz2[jk];
                rk=sqrt(k2);
                sinx=kx[ik]/rk;
                cosx=kz[jk]/rk;

                adx[ik][jk]=sinx;
                adz[ik][jk]=cosx;

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

                apx[ik][jk]=ve[0][0];
                apz[ik][jk]=ve[0][1];

                if(i!=0)
		  apvx[ik][jk]=ve[0][0]/sinx;
                else
		  apvx[ik][jk]=1.0;
			
                if(j!=0)
		  apvz[ik][jk]=ve[0][1]/cosx;
                else
		  apvz[ik][jk]=1.0;

		adx[ik][jk] *=taper[ik][jk];
		adz[ik][jk] *=taper[jk][ik];
		apx[ik][jk] *=taper[ik][jk];
		apz[ik][jk] *=taper[jk][ik];
		//apvx[ik][jk] *=taper[ik][jk];
		//apvz[ik][jk] *=taper[jk][ik];

          } /* j loop */
      } /*i loop */

     int nkx, nkz;
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
