/*************************************************************************
 * calculate finite-difference coefficents for wave propagation
 *
 * *************************************************************************/
/*
   Copyright (C) 2012 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng, Wei Kang and Tengfei Wang
     
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
#include "_fd.h"
#include "_cjb.h"

int fac(int n)
/*< fac for FD coefficent calculation >*/
{
        if(n==0)
                return 1;
        int s=1;
        while(n>0)
        {
                s=s*n;
                n--;
        }
        return s;
} 

void calc_coeff_1d(float* x)
/*< mix-th order coefficients for 1st-order spatial derivatives >*/
{
        int i;
        for(i=-_mix;i<=_mix;i++)
        {
            if(i==0)
                x[i+_mix]=0;
            else
                x[i+_mix]=2*fac(_mix)*fac(_mix)*pow(-1,i+1)/(i*fac(i+_mix)*fac(_mix-i));
        }

}

void coeff1d(float* x,float delta)
/*< m-th order coefficients for 1st-order spatial derivatives >*/
{
        int mm=2*_m+1;
        int i,j,k=0,s=0;
        float max,t,q,sum=0;

        float** A=sf_floatalloc2(mm+1,mm);

        for(i=0;i<mm;i++)
        {
                A[i][_m]=0;
                A[i][mm]=0;
        }
        A[0][_m]=1;
        A[1][mm]=1;

        for(i=0;i<mm;i++)
                for(j=-_m;j<=_m;j++)
                {
                        if(j==0)
                                continue;
                        A[i][j+_m]=pow(j*delta,i)/fac(i);
                }

        while(k<mm)
        {
                max=-9999;
                for(i=k;i<mm;i++)
                        if(A[i][k]>max)
                        {
                                max=A[i][k];
                                s=i;
                        }
                for(j=k;j<=mm;j++)
                {
                        t=A[k][j];
                        A[k][j]=A[s][j];
                        A[s][j]=t;
                }
                for(i=k+1;i<mm;i++)
                {
                        q=A[i][k]/A[k][k];
                        for(j=k;j<=mm;j++)
                                A[i][j]=A[i][j]-q*A[k][j];
                }
                k++;
        }
        for(i=mm-1;i>=0;i--)
        {
                sum=0;
                for(j=i+1;j<mm;j++)
                        sum+=A[i][j]*x[j];
                x[i]=(A[i][mm]-sum)/A[i][i];
        }
      free(*A);
}

void coeff2d(float* x,float delta)    
/*< m-th order coefficients for 2nd-order spatial derivatives, e.g. displacement-time equation >*/
{
        int mm=_m+1;       
        int i,j,k=0,s=0;
        float max,t,q,sum=0;

        float** A=sf_floatalloc2(mm+1,mm);

        for(i=0;i<mm;i++)   
        {
                A[i][mm-1]=0;
                A[i][mm]=0;
        }
        A[0][mm-1]=1;
        A[1][mm]=1;
        
        for(i=0;i<mm;i++)
                for(j=0;j<mm-1;j++)
                {
                        A[i][j]=2*pow((float)(mm-1-j)*delta,2*i)/fac(2*i);
                }
	while(k<mm)
        {
                max=-9999;
                for(i=k;i<mm;i++)
                        if(A[i][k]>max)
                        {
                                max=A[i][k];
                                s=i;
                        }                  
                for(j=k;j<=mm;j++)
                {
                        t=A[k][j];
                        A[k][j]=A[s][j];
                        A[s][j]=t;
                }                            
                for(i=k+1;i<mm;i++)
                {
                        q=A[i][k]/A[k][k];
                        for(j=k;j<=mm;j++)
                                A[i][j]=A[i][j]-q*A[k][j];
                }
                k++;
        }

        for(i=mm-1;i>=0;i--)
        {
                sum=0;
                for(j=i+1;j<mm;j++)
                        sum+=A[i][j]*x[j];
                x[i]=(A[i][mm]-sum)/A[i][i];   
                x[2*mm-2-i]=x[i];
        }

      free(*A);
}

void coeff1dmix(float* x, float delta)
/*< m-th order coefficients for 1st-order mixed spatial derivatives >*/
{
        int i;

        for(i=-_m;i<=_m;i++)
        {
                if(i==0)
                        x[i+_m]=0;
                else
                       x[i+_m]=2*fac(_m)*fac(_m)*pow(-1,i+1)/(i*fac(i+_m)*fac(_m-i));
        }
}

