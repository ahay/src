/* Functions for Chebyshev spectral method */
/*
  Copyright (C) 2008 New York University
  
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

void chebcoef(const float *b, float *coef, int n) 
/*< coefficient >*/
{
  int i,j;
  float sum;
  for(j=0;j<n;j++) {
    sum=0.0;
    for(i=0;i<n;i++) sum+=b[i]*cosf(SF_PI*j*(n-i-0.5)/n);
    coef[j]=sum*2.0/n;
  }
}

float chebeval(const float *coef, float x, int n) 
/*< evaluation >*/
{
    float *d, d0;
    int i;
    
    d = sf_floatalloc(n+2);
    
    d[n+1]=0.0;
    d[n]=0.0;
    for(i=n-1;i>0;i--) d[i]=2.0*x*d[i+1]-d[i+2]+coef[i];
    d0=x*d[1]-d[2]+0.5*coef[0];

    free(d);

    return d0;
}

void chebder(float x1, float x2, const float *coef, float *cder, int ncol) 
/*< derivative >*/
{
  int j;
  float con;

  cder[ncol-1]=0.0;
  cder[ncol-2]=2*(ncol-1)*coef[ncol-1];
  for( j=ncol-3; j>=0; j-- ) 
    cder[j]=cder[j+2]+2*(j+1)*coef[j+1];
  con=2.0/(x2-x1);
  for(j=0;j<ncol;j++) cder[j]*=con;
}

void spline(const float *x, const float *b, int n, float yp1,float ypn, float *b2) 
/*< spline coefficients >*/
{
  int i,k;
  float aux,qn,sig,un,*u;

  u= sf_floatalloc(n);
  if( yp1>0.99e30 ) b2[0]=u[0]=0.0;
  else {
    b2[0]=-0.5;
    u[0]=(3.0/(x[1]-x[0]))*((b[1]-b[0])/(x[1]-x[0])-yp1);
  }
  for( i=1; i<n-1; i++ ) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    aux=sig*b2[i-1]+2.0;
    b2[i]=(sig-1.0)/aux;
    u[i]=(b[i+1]-b[i])/(x[i+1]-x[i])-(b[i]-b[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/aux;
  }
  if( ypn>0.99e30 ) qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(b[n-1]-b[n-2])/(x[n-1]-x[n-2]));
  }
  b2[n-1]=(un-qn*u[n-2])/(qn*b2[n-2]+1.0);
  for( k=n-2; k>=0; k-- ) b2[k]=b2[k]*b2[k+1]+u[k];
  free(u);
}

float splineval(const float *x, const float *b, const float *b2, int n,float xstar) 
/*< spline evaluation >*/
{
  float h,x1,x2;
  int klo=0,khi=n-1,k;
  while(khi-klo>1) {
    k=(khi+klo)>>1;
    if(x[k]>xstar) khi=k;
    else klo=k;
  }
  h=x[khi]-x[klo];
  x1=(x[khi]-xstar)/h;
  x2=(xstar-x[klo])/h;
  return x1*b[klo]+x2*b[khi]+((x1*x1*x1-x1)*b2[klo]+(x2*x2*x2-x2)*b2[khi])*h*h/6.0;
}

