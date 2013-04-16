/* FFT routines. */
/*
  Copyright (C) 2006 The Board of Trustees of Stanford University
  
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

#include "fft.h"

void fftshift(sf_complex *a,int n)
{
  sf_complex tmp;
  //sf_complex *b;
  //b=allocatec(n);
  int n0=(n-1)/2;
  int i;
//  if (n%2==0){
    n0=n/2;
    for (i=0;i<n0;i++){
      tmp=a[i];
      a[i]=a[n0+i];
      a[n0+i]=tmp;
    }
//  }
/*  else{
    for(i=0;i<n;i++) b[i]=a[i];
    n0=n/2;
    for(i=0;i<=n0;i++)a[i]=b[i+n0];
    for(i=1;i<=n0;i++) a[n0+i]=b[i-1];
  }
*/
  //free(b);
}




void cwpfft(void *fz, int n, int isign)
{ 
  register int i; 
  register float  scale;
  float  *fp  = (float *)fz; 
  sf_complex *cz = (sf_complex *)fz;
  
  //printf("in ccc1 max0=%f,%d\n",maxvalc(cz,n),n);
  pfacc(isign, n, cz);
  //printf("in ccc2 max0=%f\n",maxvalc(cz,n));
  if( isign < 0 ) {
    scale = 1.0/(float)(n); 
    for(i=0; i < 2*n; i++) 
    {
       *fp *= scale; fp++;
    }
  }
}

void  cwpfft1d(sf_complex *a,int n,int inv)
{
  if(inv==-1){
    fftshift(a,n);
    cwpfft(a,n,-1);
  }
  else{
  //printf("in before fft max0=%f\n",maxvalc(a,n));
    cwpfft(a,n,1);
  //printf("in end fft 1 max0=%f\n",maxvalc(a,n));
    fftshift(a,n);
  //printf("in end fft 2 max0=%f\n",maxvalc(a,n));
  }
}      

int fftn(int nmin)
/*< next size >*/
{
   return npfao(nmin,nmin+50);
}
