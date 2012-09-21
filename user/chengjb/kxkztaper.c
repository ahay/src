/*************************************************************************
 * Calculate kx and kz arrays for wavenumber
 * calculate tapers on kx and kx axis
 *
 *    Copyright: Tongji University (Jiubing Cheng)
 *    2012.3.6
 *    Modified on Madagascar at 2012.8.1
 *************************************************************************/
#include "_cjb.h"

/* The tapering operators is borrowed from J. Yan */

/* Joe's taper */
#define TAPER(k) (0.5*(1+cosf(k))) 

/*#define filter(k) sin(k)/k*/
//#define TAPER(k) (k!=0 ? (k<PI ? filter(k) :0) :1)

#define filter(k) 8./5.  *sin(  k)/k -  \
                  2./5.  *sin(2*k)/k +  \
                  8./105.*sin(3*k)/k -  \
                  1./140.*sin(4*k)/k 

/* Sinc 2nd order */
#define TAPER2(k) (k!=0 ? sin(  k)/k : 1)
/*#define TAPER(k) (sin(  k)/k )*/
/* Sinc 4th order */
#define TAPER4(k) (k!=0 ?              \
                  4./3.  *sin(  k)/k - \
                  1./6.  *sin(2*k)/k : 1)
/* Sinc 6th order */
#define TAPER6(k) (k!=0 ?              \
                  3./2.  *sin(  k)/k - \
                  3./10. *sin(2*k)/k + \
                  1./60. *sin(3*k)/k : 1)
/* Sinc 8th order */
#define TAPER8(k) (k!=0 ?                    \
                  8./5.  *sin(  k)/(k) -     \
                  2./5.  *sin(2*k)/(k) +     \
                  8./105.*sin(3*k)/(k) -     \
                  1./140.*sin(4*k)/(k) : 1)

/* Dave's 2nd order */
#define TAPERD(k1,k2) (k1!=0 ? sinf(k1/2)*cosf(k2/2)/k1 :  2*cosf(k2/2) )

#define TAPERG(kmag,sig) exp(-(kmag*kmag)/2/sig/sig)/2.

#define TAPERS(k) (k!=0 ? sin(  k)/k : 1)

#define TAPERK(k1,k2) (k1*k1+k2*k2 != 0 ?                               \
                       (sqrt(k1*k1+k2*k2)<PI  ? filter(k1) :0.)      \
                       :1)

/*Butterworth 3rd order low pass filter*/
#define filterB(k)       (1+cos(k))*(1+cos(k))*(1+cos(k))/(5+3*cos(2*k))
#define filterBreal(k)   4*pow(cos(k/2),4) * (-1+3*cos(k))         /(5+3*cos(2*k))
#define filterBimag(k)   4*pow(cos(k/2),3) * ( 1+3*cos(k))*sin(k/2)/(5+3*cos(2*k))

#define TAPERB(kmag)  (kmag<PI  ? filterB(kmag)  :0 )
#define TAPERBR(kmag) (kmag<PI  ? filterBreal(kmag)  :0 )
#define TAPERBI(kmag) (kmag<PI  ? filterBimag(kmag)  :0 ) 

void kxkztaper(float *kx,float *kz, float *kkx,float *kkz, float *kx2, float *kz2, float **taper,
               int nkx, int nkz, int hnkx, int hnkz, float dkx, float dkz, float kxmax, float kzmax,
               char *tapertype)
/*< kxkztaper: define kx and kz array for kx-axis and kz-aixs wavenumber, and also calculate tapers on kx and kx axis >*/
{
        int   i, j, ik, jk;
        float rkx, rkz, k2;
        float *taperx, *taperz;
        int order=8;

         for( i=-hnkx; i<=hnkx ; i++ )
         {
             ik=i+hnkx;
             for( j=-hnkz; j<=hnkz ; j++)
             {
               jk=j+hnkz;
               taper[ik][jk] = 1.0;
             }
         }
       
         if(tapertype[0]=='C'||tapertype[0]=='S'||tapertype[0]=='T')
         {
           taperx=calloc(sizeof(float), nkx);
           taperz=calloc(sizeof(float), nkz);
         }
     
         for( i=-hnkx; i<=hnkx ; i++ )
         {
           ik=i+hnkx;
           kx[ik]=i*dkx;
           kx2[ik]=kx[ik]*kx[ik];
           kkx[ik] = 2*PI*i/nkx;
           rkx=kkx[ik];
           if(tapertype[0]=='S')// Jia Yan's Sinc taper
           {
              switch(order){
                case 8: taperx[ik]=TAPER8(rkx); break;
                case 6: taperx[ik]=TAPER6(rkx); break;
                case 4: taperx[ik]=TAPER4(rkx); break;
                case 2: taperx[ik]=TAPER2(rkx);
              }
           }else if(tapertype[0]=='C'||tapertype[0]=='T')//Joe Dellinger's Cosine taper
           {
              /* 0.5*(1.0+cos(PI*kx[ik]/kx_nyquist)) */
              taperx[ik]=TAPER(rkx);
           }
        }
        for( j=-hnkz; j<=hnkz ; j++)
        {
           jk=j+hnkz;
           kz[jk]=j*dkz;
           kz2[jk]=kz[jk]*kz[jk];
           kkz[jk] = 2*PI*j/nkz;
           rkz=kkz[jk];
           if(tapertype[0]=='S') // Jia Yan's Sinc taper
           {

              switch(order){
                case 8: taperz[jk]=TAPER8(rkz); break;
                case 6: taperz[jk]=TAPER6(rkz); break;
                case 4: taperz[jk]=TAPER4(rkz); break;
                case 2: taperz[jk]=TAPER2(rkz);
              }
           }else if(tapertype[0]=='C'||tapertype[0]=='T')//Joe Dellinger's Cosine taper
           {
              /* 0.5*(1.0+cos(PI*kx[ik]/kx_nyquist)) */
              taperz[jk]=TAPER(rkz);
           }
        }

        /******************* calculate 2D taper *************************/
        if(tapertype[0]=='D')// Dale Hale's taper
        {
          float sig=1.5;          

          for( i=-hnkx; i<=hnkx ; i++ )
          {
             ik=i+hnkx;
             for( j=-hnkz; j<=hnkz ; j++)
             {
               jk=j+hnkz;
               k2=kx2[ik]+kz2[jk];
               taper[ik][jk] = TAPERG(k2, sig);
               //taper[ik][jk] = TAPERD(kx[ik], kz[jk]);
             }
           }
        }else if(tapertype[0]=='K')// 
        {
          for( i=-hnkx; i<=hnkx ; i++ )
          {
             ik=i+hnkx;
             for( j=-hnkz; j<=hnkz ; j++)
             {
               jk=j+hnkz;
               taper[ik][jk] = TAPERK(kx[ik], kz[jk]);
             }
           }
        }else if(tapertype[0]=='T') //CJB' Taper for TTI
        {
          for( i=-hnkx; i<=hnkx ; i++ )
          {
             ik=i+hnkx;
             for( j=-hnkz; j<=hnkz ; j++)
             {
               jk=j+hnkz;
               taper[ik][jk] = taperx[ik]*taperz[jk];
             }
           }
        }else //Joe Dellinger's Cosine taper or Jia Yan's Sinc Taper
        {
          for( i=-hnkx; i<=hnkx ; i++ )
          {
             ik=i+hnkx;
             for( j=-hnkz; j<=hnkz ; j++)
             {
               jk=j+hnkz;
               taper[ik][jk] = taperx[ik];
             }
           }
        }

        if(tapertype[0]=='C'||tapertype[0]=='S'||tapertype[0]=='T')
        {
           free(taperx);
           free(taperz);
        }
}
