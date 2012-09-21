/*************************************************************************
 * Calculate projection deviation operators in wavenumber domain 
 * 1) define wave-vector and devergence operators;
 * 2) solve kelvin-Christoffel equation for the eigenvalue & eigeinvectors
 * 3) calculate projection deviation due to difference between 
 *    wave-vector and P-wave's polarization vector
 *
 *    Copyright: Tongji University (Jiubing Cheng)
 *    2012.3.2
 *************************************************************************/
#include <rsf.h>

#include "_cjb.h"
#include "smth2d.h"
#include "engein2dvti.h"

void devttip(float **apvx,float **apvz,
             float *kx,float *kz, float *kkx,float *kkz, float **taper,          
             int hnkx, int hnkz,
             float vp2, float vs2, float ep2, float de2, float the, int itaper, int ismth)
/*< devttip: projection deviation operators for P-wave in TTI media >*/
{
        int   i, j, ii, jj, ik, jk, num;
        float kx2, kz2, k2, rk, sinx, cosx;
        float coss, sins;
        float sum, kxx, kzz;

        float ve[2][2],va[2];

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

                //engein2dvti1(ve, va, sinx, cosx, vp2, vs2, ep2, de2, f);
                
                engein2dvti2(ve, va, sinx, cosx, vp2, vs2, ep2, de2);
                
                //engein2dvti3(ve, va, sinx, cosx, vp2, vs2, ep2, de2);

                if(i!=0)
		  apvx[ik][jk]=ve[0][0]/sinx;
                else
		  apvx[ik][jk]=1.0;
			
                if(j!=0)
		  apvz[ik][jk]=ve[0][1]/cosx;
                else
		  apvz[ik][jk]=1.0;

		apvx[ik][jk] *=taper[ik][jk];
		apvz[ik][jk] *=taper[jk][ik];

          } /* j loop */
      } /*i loop */

     int nkx, nkz;
     nkx=2*hnkx+1;
     nkz=2*hnkz+1;

     /* interpolating */
     for( i=0; i<nkx; i++ )
     {
          apvz[i][hnkz]=(apvz[i][hnkz+1] + apvz[i][hnkz-1])/2.0;
          apvx[i][hnkz]=(apvx[i][hnkz+1] + apvx[i][hnkz-1])/2.0;
     }
     for( j=0; j<nkz ; j++)
     {
          apvx[hnkx][j]=(apvx[hnkx+1][j] + apvx[hnkx-1][j])/2.0;
          apvz[hnkx][j]=(apvz[hnkx+1][j] + apvz[hnkx-1][j])/2.0;
     }

     /* interpolating & smoothing for PI/4 */
     if(fabs(fabs(the)-PI/4.0)<0.001)
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

     /***************************************************************************
     void smooth2d(float **v, int n1, int n2, float r1, float r2, float rw)
     **************************************************************************
     int n1;          number of points in x1 (fast) dimension
     int n2;          number of points in x1 (fast) dimension 
     float **v0;       array of input velocities 
     float r1;        smoothing parameter for x1 direction
     float r2;        smoothing parameter for x2 direction
     float rw;        smoothing parameter for window
     " Notes:                                                                ",
     " Larger r1 and r2 result in a smoother data. Recommended ranges of r1  ", 
     " and r2 are from 1 to 20.                                              ",
     "                                                                       ",
     " Smoothing can be implemented in a selected window. The range of 1st   ",
     " dimension for window is from win[0] to win[1]; the range of 2nd       ",
     " dimension is from win[2] to win[3].                                   ",
     "                                                                       ",
     " Smoothing the window function (i.e. blurring the edges of the window) ",
     " may be done by setting a nonzero value for rw, otherwise the edges    ",
     " of the window will be sharp.                                          ",
     "                                                                       ",
     ***************************************************************************/

     float r1, r2, rw;

     r1=4.0;
     r2=4.0;
     rw=2.0;

     if(ismth==1)
     {
       smooth2d(apvx, nkz, nkx, r1, r2, rw);
       smooth2d(apvz, nkz, nkx, r1, r2, rw);
     }
}
