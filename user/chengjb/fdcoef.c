/*************************************************************************
 * calculate finite-difference coefficents for wave propagation
 *
 * Copyright: Tongji University (Jiubing Cheng)
 * 2012.3.2
 *
 * *************************************************************************/
#include <rsf.h>
#include "_fd.h"
/* #include "_cjb.h" */
#include "alloc.h"

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

void coeff1d(float* x,float delta)
/*< m-th order coefficients for 1st-order spatial derivatives >*/
{
        int mm=2*m+1;
        int i,j,k=0,s=0;
        float max,t,q,sum=0;

        float** A=sf_floatalloc2(mm+1,mm);

        for(i=0;i<mm;i++)
        {
                A[i][m]=0;
                A[i][mm]=0;
        }
        A[0][m]=1;
        A[1][mm]=1;

        for(i=0;i<mm;i++)
                for(j=-m;j<=m;j++)
                {
                        if(j==0)
                                continue;
                        A[i][j+m]=pow(j*delta,i)/fac(i);
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
      free2float(A);
}

void coeff2d(float* x,float delta)    
/*< m-th order coefficients for 2nd-order spatial derivatives >*/
{
        int mm=m+1;       
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

      free2float(A);
}

void coeff1dmix(float* x,float delta)
/*< m-th order coefficients for 1st-order mixed spatial derivatives >*/
{
        int i;

        for(i=-m;i<=m;i++)
        {
                if(i==0)
                        x[i+m]=0;
                else
                       x[i+m]=2*fac(m)*fac(m)*pow(-1,i+1)/(i*fac(i+m)*fac(m-i));
        }
}

