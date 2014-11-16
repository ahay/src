/* Apply convolution kernel or its adjoint (without cyclic boundaries)*/
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include <math.h>
#include <float.h>

#include <rsf.h>
/*^*/

#include "convkernel.h"

static int n1,n2,n3,nl;
static int *lag1, *lag2, *lag3;
static float *filter;

void convkernel1_init(int n1_in /* data size */,
                      int nl_in /* number of filter lags */,
                      int *lag1_in /* filter lags */,
                      float* filter_in /* filter coeficients */) 
/*< initialize >*/
{
    n1 = n1_in;
    nl = nl_in;
    lag1 = lag1_in;
    filter = filter_in;
}

void convkernel2_init(int n1_in /* data size along axis 1 */,
                      int n2_in /* data size along axis 2 */,
                      int nl_in /* number of filter lags */,
                      int *lag1_in /* filter lags along axis 1*/,
                      int *lag2_in /* filter lags along axis 2*/,
                      float *filter_in /* filter coeficients */) 
/*< initialize2 >*/
{
    n1 = n1_in;
    n2 = n2_in;
    nl = nl_in;
    lag1 = lag1_in;
    lag2 = lag2_in;
    filter = filter_in;
}


void convkernel3_init(int n1_in /* data size along axis 1 */,
                      int n2_in /* data size along axis 2 */,
                      int n3_in /* data size along axis 3 */,
                      int nl_in /* number of filter lags */,
                      int *lag1_in /* filter lags along axis 1*/,
                      int *lag2_in /* filter lags along axis 2*/,
                      int *lag3_in /* filter lags along axis 3*/,
                      float *filter_in /* filter coeficients */) 
/*< initialize3 >*/
{
    n1 = n1_in;
    n2 = n2_in;
    n3 = n3_in;
    nl = nl_in;
    lag1 = lag1_in;
    lag2 = lag2_in;
    lag3 = lag3_in;
    filter = filter_in;
}




void convkernel1_apply(
        float *xx,
        float *yy,
        bool adj)
/*< apply kernel1 >*/ 
{
  int i1; 
  int ilag;
  int i_1;

  for (i1=0; i1<n1; ++i1){
    if(adj){
      xx[i1] = 0.0;
    }else{
      yy[i1] = 0.0;
    }
  }



  for (i1=0; i1<n1; ++i1){
    for (ilag=0; ilag<nl; ++ilag){
      i_1 = i1 - lag1[ilag];
      if(i_1 <=0 || i_1 >=n1 || i1 <= 0 || i1>=n1) continue;

      if(adj){
        xx[i_1] += filter[ilag]*yy[i1 ];
      }else{
        yy[i1 ] += filter[ilag]*xx[i_1];
      }
    }
  }
}



void convkernel2_apply(
        float **xx,
        float **yy,
        bool adj)
/*< apply kernel2 >*/ 
{
  int i1,i2; 
  int ilag;
  int i_1,i_2;

  for (i2=0; i2<n2; ++i2){
    for (i1=0; i1<n1; ++i1){
      if(adj){
        xx[i2][i1] = 0.0;
      }else{
        yy[i2][i1] = 0.0;
      }
    }
  }



  for (ilag=0; ilag<nl; ++ilag){
    for (i2=0; i2<n2; ++i2){
      i_2 = i2 - lag2[ilag];
      if(i_2 <=0 || i_2 >=n2 || i2 <= 0 || i2>=n2) continue;

      for (i1=0; i1<n1; ++i1){
        i_1 = i1 - lag1[ilag];
        if(i_1 <=0 || i_1 >=n1 || i1 <= 0 || i1>=n1) continue;
  
        if(adj){
          xx[i_2][i_1] += filter[ilag]*yy[i2 ][i1 ];
        }else{
          yy[i2 ][i1 ] += filter[ilag]*xx[i_2][i_1];
        }
      }
    }
  }
}

void convkernel3_apply(
        float ***xx,
        float ***yy,
        bool adj)
/*< apply kernel3 >*/ 
{
  int i1,i2,i3; 
  int ilag;
  int i_1,i_2,i_3;

  for (i3=0; i3<n3; ++i3){
    for (i2=0; i2<n2; ++i2){
      for (i1=0; i1<n1; ++i1){
        if(adj){
          xx[i3][i2][i1] = 0.0;
        }else{
          yy[i3][i2][i1] = 0.0;
        }
      }
    }
  }

  for (ilag=0; ilag<nl; ++ilag){
    for (i3=0; i3<n3; ++i3){
      i_3 = i3 - lag3[ilag];
      if(i_3 <=0 || i_3 >=n3 || i3 <= 0 || i3>=n3) continue;

      for (i2=0; i2<n2; ++i2){
        i_2 = i2 - lag2[ilag];
        if(i_2 <=0 || i_2 >=n2 || i2 <= 0 || i2>=n2) continue;

        for (i1=0; i1<n1; ++i1){
          i_1 = i1 - lag1[ilag];
          if(i_1 <=0 || i_1 >=n1 || i1 <= 0 || i1>=n1) continue;
    
          if(adj){
            xx[i_3][i_2][i_1] += filter[ilag]*yy[i3 ][i2 ][i1 ];
          }else{
            yy[i3 ][i2 ][i1 ] += filter[ilag]*xx[i_3][i_2][i_1];
          }
        }
      }
    }
  }
}


