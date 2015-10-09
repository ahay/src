/*************************************************************************
 * Calculate kx and kz arrays for wavenumber
 * calculate tapers on kx and kx axis
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

/* The tapering operators is borrowed from J. Yan */

/* Joe's taper */
#define TAPER(k) (0.5*(1+cosf(k))) 

/*#define filter(k) sin(k)/k*/
/*#define TAPER(k) (k!=0 ? (k<PI ? filter(k) :0) :1)*/

#define filter(k) 8./5.  *sin(  k)/k -			\
				  2./5.  *sin(2*k)/k +  \
				  8./105.*sin(3*k)/k -  \
				  1./140.*sin(4*k)/k 

/* Sinc 2nd order */
#define TAPER2(k) (k!=0 ? sin(  k)/k : 1)
/*#define TAPER(k) (sin(  k)/k )*/
/* Sinc 4th order */
#define TAPER4(k) (k!=0 ?			\
		   4./3.  *sin(  k)/k -		\
		   1./6.  *sin(2*k)/k : 1)
/* Sinc 6th order */
#define TAPER6(k) (k!=0 ?			\
		   3./2.  *sin(  k)/k -		\
		   3./10. *sin(2*k)/k +		\
		   1./60. *sin(3*k)/k : 1)
/* Sinc 8th order */
#define TAPER8(k) (k!=0 ?			\
		   8./5.  *sin(  k)/(k) -	\
		   2./5.  *sin(2*k)/(k) +	\
		   8./105.*sin(3*k)/(k) -	\
		   1./140.*sin(4*k)/(k) : 1)

/* Dave's 2nd order */
#define TAPERD(k1,k2) (k1!=0 ? sinf(k1/2)*cosf(k2/2)/k1 :  2*cosf(k2/2) )

#define TAPERG(kmag,sig) exp(-(kmag*kmag)/2/sig/sig)/2.

#define TAPERS(k) (k!=0 ? sin(  k)/k : 1)

#define TAPERK(k1,k2) (k1*k1+k2*k2 != 0 ?                               \
                       (sqrt(k1*k1+k2*k2)<SF_PI  ? filter(k1) :0.)      \
                       :1)

/*Butterworth 3rd order low pass filter*/
#define filterB(k)       (1+cos(k))*(1+cos(k))*(1+cos(k))/(5+3*cos(2*k))
#define filterBreal(k)   4*pow(cos(k/2),4) * (-1+3*cos(k))         /(5+3*cos(2*k))
#define filterBimag(k)   4*pow(cos(k/2),3) * ( 1+3*cos(k))*sin(k/2)/(5+3*cos(2*k))

#define TAPERB(kmag)  (kmag<SF_PI  ? filterB(kmag)  :0 )
#define TAPERBR(kmag) (kmag<SF_PI  ? filterBreal(kmag)  :0 )
#define TAPERBI(kmag) (kmag<SF_PI  ? filterBimag(kmag)  :0 ) 

void ikx(int *ijkx, int nkx)
/*< ikx: re-define the index for kx axis after FFT >*/
{
    int ix;
 
    if(nkx%2==0){
	for(ix=0;ix<nkx;ix++)
	{
	    if(ix<nkx/2)
		ijkx[ix]=nkx/2+ix;
	    else
		ijkx[ix]=ix-nkx/2;
	}
    }else{
	for(ix=0;ix<nkx;ix++)
	{
	    if(ix<=nkx/2)
		ijkx[ix]=nkx/2+ix;
	    else
		ijkx[ix]=ix-nkx/2-1;
	}
    }
}

void ikxikz(int *ijkx, int *ijkz, int nkx, int nkz)
/*< ikxikz: re-define the index for kx and kz axis after FFT >*/
{
    int ikx, ikz;
 
    if(nkx%2==0){
	for(ikx=0;ikx<nkx;ikx++)
	{
	    if(ikx<nkx/2)
		ijkx[ikx]=nkx/2+ikx;
	    else
		ijkx[ikx]=ikx-nkx/2;
	}
    }else{
	for(ikx=0;ikx<nkx;ikx++)
	{
	    if(ikx<=nkx/2)
		ijkx[ikx]=nkx/2+ikx;
	    else
		ijkx[ikx]=ikx-nkx/2-1;
	}
    }

    if(nkz%2==0){
	for(ikz=0;ikz<nkz;ikz++)
	{
	    if(ikz<nkz/2)
		ijkz[ikz]=nkz/2+ikz;
	    else
		ijkz[ikz]=ikz-nkz/2;
	}
    }else{
	for(ikz=0;ikz<nkz;ikz++)
	{
	    if(ikz<=nkz/2)
		ijkz[ikz]=nkz/2+ikz;
	    else
		ijkz[ikz]=ikz-nkz/2-1;
	}
    }
}

void ikxikyikz(int *ijkx, int *ijky, int *ijkz, int nkx, int nky, int nkz)
/*< ikxikyikz: re-define the index for kx and kz axis after FFT >*/
{
    int ikx, iky, ikz;
 
    if(nky%2==0){
	for(iky=0;iky<nky;iky++)
	{
	    if(iky<nky/2)
		ijky[iky]=nky/2+iky;
	    else
		ijky[iky]=iky-nky/2;
	}
    }else{
	for(iky=0;iky<nky;iky++)
	{
	    if(iky<=nky/2)
		ijky[iky]=nky/2+iky;
	    else
		ijky[iky]=iky-nky/2-1;
	}
    }

    if(nkx%2==0){
	for(ikx=0;ikx<nkx;ikx++)
	{
	    if(ikx<nkx/2)
		ijkx[ikx]=nkx/2+ikx;
	    else
		ijkx[ikx]=ikx-nkx/2;
	}
    }else{
	for(ikx=0;ikx<nkx;ikx++)
	{
	    if(ikx<=nkx/2)
		ijkx[ikx]=nkx/2+ikx;
	    else
		ijkx[ikx]=ikx-nkx/2-1;
	}
    }

    if(nkz%2==0){
	for(ikz=0;ikz<nkz;ikz++)
	{
	    if(ikz<nkz/2)
		ijkz[ikz]=nkz/2+ikz;
	    else
		ijkz[ikz]=ikz-nkz/2;
	}
    }else{
	for(ikz=0;ikz<nkz;ikz++)
	{
	    if(ikz<=nkz/2)
		ijkz[ikz]=nkz/2+ikz;
	    else
		ijkz[ikz]=ikz-nkz/2-1;
	}
    }
}

void kxkztaper(float *kx,float *kz, float *kkx,float *kkz, float *kx2, float *kz2, float **taper,
               int nkx, int nkz, int hnkx, int hnkz, float dkx, float dkz, float kxmax, float kzmax,
               char *tapertype)
/*< kxkztaper: define kx and kz array for kx-axis and kz-aixs wavenumber, and also calculate tapers on kx and kx axis >*/
{
    int   i, j, ik, jk;
    float rkx, rkz, k2;
    float *taperx, *taperz;
    int order=8;
    float sig=1.5;          

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
	taperx=sf_floatalloc(nkx);
	taperz=sf_floatalloc(nkz);
    } else {
	taperx=NULL;
	taperz=NULL;
    }
     
    for( i=-hnkx; i<=hnkx ; i++ )
    {
	ik=i+hnkx;
	kx[ik]=i*dkx;
	kx2[ik]=kx[ik]*kx[ik];
	kkx[ik] = 2*SF_PI*i/nkx;
	rkx=kkx[ik];
	if(tapertype[0]=='S')/* Jia Yan's Sinc taper */
	{
	    switch(order){
                case 8: taperx[ik]=TAPER8(rkx); break;
                case 6: taperx[ik]=TAPER6(rkx); break;
                case 4: taperx[ik]=TAPER4(rkx); break;
                case 2: taperx[ik]=TAPER2(rkx);
	    }
	}else if(tapertype[0]=='C'||tapertype[0]=='T')/*Joe Dellinger's Cosine taper */
	{
	    /* 0.5*(1.0+cos(SF_PI*kx[ik]/kx_nyquist)) */
	    taperx[ik]=TAPER(rkx);
	}
    }
    for( j=-hnkz; j<=hnkz ; j++)
    {
	jk=j+hnkz;
	kz[jk]=j*dkz;
	kz2[jk]=kz[jk]*kz[jk];
	kkz[jk] = 2*SF_PI*j/nkz;
	rkz=kkz[jk];
	if(tapertype[0]=='S') /* Jia Yan's Sinc taper */
	{

	    switch(order){
                case 8: taperz[jk]=TAPER8(rkz); break;
                case 6: taperz[jk]=TAPER6(rkz); break;
                case 4: taperz[jk]=TAPER4(rkz); break;
                case 2: taperz[jk]=TAPER2(rkz);
	    }
	}else if(tapertype[0]=='C'||tapertype[0]=='T')/*Joe Dellinger's Cosine taper */
	{
	    /* 0.5*(1.0+cos(SF_PI*kx[ik]/kx_nyquist)) */
	    taperz[jk]=TAPER(rkz);
	}
    }

    /******************* calculate 2D taper *************************/
    if(tapertype[0]=='D')/* Dale Hale's taper */
    {
	for( i=-hnkx; i<=hnkx ; i++ )
	{
	    ik=i+hnkx;
	    for( j=-hnkz; j<=hnkz ; j++)
	    {
		jk=j+hnkz;
		k2=kx2[ik]+kz2[jk];
		taper[ik][jk] = TAPERG(k2, sig);
		/* taper[ik][jk] = TAPERD(kx[ik], kz[jk]); */
	    }
	}
    }else if(tapertype[0]=='K') 
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
    }else if(tapertype[0]=='T') /* CJB' Taper for TTI */
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
    }else /* Joe Dellinger's Cosine taper or Jia Yan's Sinc Taper */
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

void kykxkztaper(float *ky, float *kx,float *kz, float *kky, float *kkx,float *kkz,
                 float *ky2, float *kx2, float *kz2, float ***taper,
                 int nky, int nkx, int nkz, int hnky, int hnkx, int hnkz,
                 float dky, float dkx, float dkz, float kymax, float kxmax, float kzmax)
/*< kykxkztaper: define ky, kx and kz array for ky-axis, kx-axis and kz-aixs wavenumber,
 * and also calculate tapers on ky, kx and kx axis >*/
{
    int   i, j, k, ik, jk, kk;
    float rky, rkx, rkz;
    float *taperx, *tapery, *taperz;

    for( k=-hnky; k<=hnky ; k++ )
    {
	kk=k+hnky;
	for( i=-hnkx; i<=hnkx ; i++ )
	{
	    ik=i+hnkx;
	    for( j=-hnkz; j<=hnkz ; j++)
	    {
		jk=j+hnkz;
		taper[kk][ik][jk] = 1.0;
	    }
	}
    }
       
    tapery=sf_floatalloc(nky);
    taperx=sf_floatalloc(nkx);
    taperz=sf_floatalloc(nkz);
     
    for( i=-hnkx; i<=hnkx ; i++ )
    {
	ik=i+hnkx;
	kx[ik]=i*dkx;
	kx2[ik]=kx[ik]*kx[ik];
	kkx[ik] = 2*SF_PI*i/nkx;
	rkx=kkx[ik];
	/* 0.5*(1.0+cos(SF_PI*kx[ik]/kx_nyquist)) */
	taperx[ik]=TAPER(rkx);
    }
    for( j=-hnkz; j<=hnkz ; j++)
    {
	jk=j+hnkz;
	kz[jk]=j*dkz;
	kz2[jk]=kz[jk]*kz[jk];
	kkz[jk] = 2*SF_PI*j/nkz;
	rkz=kkz[jk];
	/* 0.5*(1.0+cos(SF_PI*kx[ik]/kx_nyquist)) */
	taperz[jk]=TAPER(rkz);
    }
    for( k=-hnky; k<=hnky ; k++)
    {
	kk=k+hnky;
	ky[kk]=k*dky;
	ky2[kk]=ky[kk]*ky[kk];
	kky[kk] = 2*SF_PI*k/nky;
	rky=kky[kk];
	/* 0.5*(1.0+cos(SF_PI*kx[ik]/kx_nyquist)) */
	tapery[kk]=TAPER(rky);
    }

    /******************* calculate 3D taper *************************/
    for( k=-hnky; k<=hnky ; k++ )
    {
	kk=k+hnky;
	for( i=-hnkx; i<=hnkx ; i++ )
	{
	    ik=i+hnkx;
	    for( j=-hnkz; j<=hnkz ; j++)
	    {
		jk=j+hnkz;
		taper[kk][ik][jk] = taperx[ik]*tapery[kk]*taperz[jk];
	    }
	}
    }
    free(tapery);
    free(taperx);
    free(taperz);
}
