/*************************************************************************
 * Calculate polarization operators in wavenumber domain 
 * 1) solve kelvin-Christoffel equation for the eigenvalue & eigeinvectors
 * 2) calculate P-wave's polarization operators 
 * 3) calculate SV-wave's polarization operators 
 *
 *    Copyright: Tongji University (Jiubing Cheng)
 *    2012.3.6
 *    Modified on Madagascar at 2012.8.1
 *************************************************************************/
#include <rsf.h>

#include "_cjb.h"
#include "smth2d.h"
#include "engein2dvti.h"

void polvtipsv(float **apx,float **apz, float **apxs,float **apzs,
             float *kx,float *kz, float *kkx,float *kkz, float *kx2, float *kz2,
             float **taper, int hnkx, int hnkz,
             float vp2, float vs2, float ep2, float de2, int itaper)
/*< polvtipsv: P- and SV-wave polarization operators in VTI media >*/
{
        int   i, j, ik, jk, ismth=0;
        float k2, rk, sinx, cosx;

        //float e1=0, e2=0, e12;
        float r1, r2, rw;
        float ve[2][2], va[2];
 
        //float f = 1.0-vs2/vp2;

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

                k2=kx2[ik]+kz2[jk];
                rk=sqrt(k2);
                sinx=kx[ik]/rk;
                cosx=kz[jk]/rk;

                //engein2dvti1(ve, va, sinx, cosx, vp2, vs2, ep2, de2, f);
                engein2dvti2(ve, va, sinx, cosx, vp2, vs2, ep2, de2);
                //engein2dvti3(ve, va, sinx, cosx, vp2, vs2, ep2, de2);

                apx[ik][jk]=ve[0][0];
                apz[ik][jk]=ve[0][1];
                if(ve[0][0]*sinx + ve[0][1]*cosx <=0.0)
                {
                     apx[ik][jk]= -ve[0][0];
                     apz[ik][jk]= -ve[0][1];
                }
                apx[ik][jk] *= taper[ik][jk];
                apz[ik][jk] *= taper[jk][ik];

                apxs[ik][jk]=ve[1][0];
                apzs[ik][jk]=ve[1][1];
                if(ve[1][0]*cosx - ve[1][1]*sinx <=0.0)
                {
                     apxs[ik][jk]= -ve[1][0];
                     apzs[ik][jk]= -ve[1][1];
                }
                apxs[ik][jk] *= taper[jk][ik];
                apzs[ik][jk] *= taper[ik][jk];
          } /* j loop */
      } /*i loop */

     if(ismth==1)
     {

        int nkx, nkz;
        nkx=2*hnkx+1;
        nkz=2*hnkz+1;

        r1=6.0;
        r2=6.0;
        rw=2.0;
        smooth2d(apx, nkz, nkx, r1, r2, rw);
        smooth2d(apz, nkz, nkx, r1, r2, rw);
        smooth2d(apxs, nkz, nkx, r1, r2, rw);
        smooth2d(apzs, nkz, nkx, r1, r2, rw);
     }
}
