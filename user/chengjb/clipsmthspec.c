/*************************************************************************
 * clip specified value and smoothing the spectrum
 * *************************************************************************/
#include <rsf.h>

#include "_cjb.h"

void clipsmthspec(float **a,int nx,int nz, float var)
/*<clipsmthspec: clip the values larger than "var" by smoothing >*/
{
      int i, j, ii, jj, num;
      float s;
      float **tmp=sf_floatalloc2(nz, nx);

      for( i=0; i<nx; i++ )
         for( j=0; j<nz; j++)
              tmp[i][j]=a[i][j];

      for( i=0; i<nx; i++ )
         for( j=0; j<nz; j++)
              if(fabs(tmp[i][j])>=var)
              {
                 num=0;
                 s=0.0;
                 for(ii=i-2;ii<=i+2;ii++)
                    if(ii>=0&&ii<nx)
                        for(jj=j-2;jj<=j+2;jj++)
                         if(jj>=0&&jj<nz&&i!=ii&&j!=jj&&(jj!=nz/2&&ii!=nx/2))
                          {
                              num++;
                              s += tmp[ii][jj];
                          }
                          tmp[i][j]=s/num;
               }

        for( i=0; i<nx; i++ )
           for( j=0; j<nz ; j++)
               a[i][j]=tmp[i][j];

        free(*tmp);
}
