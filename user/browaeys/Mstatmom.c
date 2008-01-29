/**** FILE: statisticlib.c                                                         ****/
/**** Statistical description of data moments                                      ****/
/****                                                                              ****/
/**** moment(data[],n,*ave,*adev,*sdev,*var,*skew,*curt)   : distribution moments  ****/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "statisticlib.h"

void moment(float data[], int n, float *ave, float *adev,
float *sdev, float *var, float *skew, float *curt)
/* Input is array data[0..n-1] and n and output are mean (ave), mean absolute deviation */
/* (adev), standard deviation (sdev), variance (var), skewness (skew) and kurtosis (curt) */
{
 int j;
 float ep, s, p;
 ep = 0.0;
 if (n <= 1) 
 {
  printf("Error in function moment, n must be at least 2 \n");
  exit(1);
 }
 s = 0.0;
 for (j=0;j<n;j++) s += data[j];
 *ave = s/n;
 *adev = (*var) = (*skew) = (*curt) = 0.0;
 for (j=0;j<n;j++)
 {
  *adev += fabs(s=data[j]-(*ave));
  ep += s;
  *var += (p=s*s);
  *skew += (p *= s);
  *curt += (p *= s);
 }
 *adev /= n;
 *var = (*var-ep*ep/n)/(n-1);
 *sdev = sqrt(*var);
 if (*var)
 {
  *skew /= (n*(*var)*(*sdev));
  *curt = (*curt)/(n*(*var)*(*var))-3.0;
 }
 else
 {
  printf(" No skew/kurtosis when variance = 0 (in function moment) \n");
 }
}
