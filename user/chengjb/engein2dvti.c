/*************************************************************************
 * Calculate eigeinvalues and eigeinvectors for 2D VTI media
 * Method 1: analytic solution based on Dellinger's expression
 * Method 2: using 2*2 matrix analytic solution
 * Method 3: using general N*N matrix Lapack solution
 *
 * Note:  Method 2 & 3 are similar, they solve same equation using different codes
 *
 *    Copyright: Tongji University (Jiubing Cheng)
 *    2012.3.2
 *************************************************************************/
#include <rsf.h>

#include "_cjb.h"
#include "_lapack.h"
#include "eigen2x2.h"

void engein2dvti1(float ve[2][2], float va[2], float sinx, float cosx, float vp2, float vs2,
                 float ep2, float de2, float f)
/*< engein2dvti1: Calculate eigeinvalues and eigeinvectors for 2D VTI media
                  using analytic solution based on Dellinger's expression >*/
{
        float sin2, cos2, cos2a, sin4;
        float d33, d11, psi2, psi, d33d11, d33d11_2, sin2cos2, sin2cos2_4, tmpa,tmpb,u1,u2,u1u2;

        sin2=sinx*sinx;
        cos2=1.0-sin2;
        cos2a=1.0-2*sin2;
        sin4=sin2*sin2;
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

void engein2dvti2(float ve[2][2], float va[2], float kx, float kz, float vp2, float vs2,
                 float ep2, float de2)
/*< engein2dvti2: Calculate eigeinvalues and eigeinvectors for 2D VTI media
                  using 2*2 matrix analytic solution>*/
{
        float u1, u2, c11, c33, c44, c13, a11, a12, a22;
        float a[2][2];

        c33=vp2;
        c44=vs2;
        c11=ep2*c33;
        //c13=sqrt(2*c33*(c33-c44)*de+(c33-c44)*(c33-c44))-c44;
        //c13=sqrt((2*c33*de+(c33-c44))*(c33-c44))-c44;
        //c13=sqrt((2*c33*de+c33-c44)*(c33-c44))-c44;
        //c13=sqrt((2*de+1.0)*c33-c44)*(c33-c44))-c44;
        c13=sqrt((de2*c33-c44)*(c33-c44))-c44;

        a11=  c11*kx*kx+
              c44*kz*kz;

        a12= (c13+c44)*kx*kz;

        a22=  c44*kx*kx+
              c33*kz*kz;

         a[0][0] = a11;
         a[0][1] = a12;
         a[1][0] = a12;
         a[1][1] = a22;

        solveSymmetric22(a, ve, va);

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

void engein2dvti3(float ve[2][2], float va[2], float sinx, float cosx, float vp2, float vs2,
                 float ep2, float de2)
/*< engein2dvti3: Calculate eigeinvalues and eigeinvectors for 2D VTI media
                  using general N*N matrix Lapack solution >*/
{
        float sin2, cos2;

	char	jobz='V';  /* for SVD */
	char	uplo='U';  /* for SVD */
	int     M=2;  /* for SVD */
	int     LDA=M;  /* for SVD */
	int     LWORK=4*M;  /* for SVD */

	int     INFO;  /* for SVD */

        float   *Chr, *w, *work;  /* Lapack SVD array */
	Chr=calloc(sizeof(float),M*M);
        w=calloc(sizeof(float),M);
	work=calloc(sizeof(float),LWORK);
        
        sin2=sinx*sinx;
        cos2=1-sin2;

        Chr[0]=ep2*vp2*sin2+vs2*cos2;
        Chr[1]=sqrt((vp2-vs2)*(de2*vp2-vs2))*sinx*cosx;
        Chr[2]=Chr[1];
        Chr[3]=vs2*sin2+vp2*cos2;

        ssyev_(&jobz, &uplo, &M, Chr, &LDA, w, work, &LWORK, &INFO);

        /* large eigeinvalue */
        va[0]=w[1];
        ve[0][0]=Chr[2];
        ve[0][1]=Chr[3];
        if(Chr[2]*sinx + Chr[3]*cosx <=0.0)
        {
             ve[0][0]= -Chr[2];
             ve[0][1]= -Chr[3];
        }
        /* small eigeinvalue */
        va[1]=w[0];
        ve[1][0]=Chr[0];
        ve[1][1]=Chr[1];
        if(Chr[0]*cosx - Chr[1]*sinx <=0.0)
        {
             ve[1][0]= -Chr[0];
             ve[1][1]= -Chr[1];
        }

        free(Chr);
        free(w);
        free(work);
}
