/**** FILE: variolib.c                                                                   ****/
/**** Calculations of spatial geostatistics                                              ****/
/****                                                                                    ****/
/**** svariogram(nbin,*number,**data,*v)                : isotropic sample variogram     ****/
/**** smadogram(nbin,*number,**data,*v)                 : isotropic madogram estimator   ****/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "variolib.h"

void svariogram(int nbin, int *number, float **data, float *v)
/* Isotropic sample variogram calculation */
/* The discrete distance between pair of locations is binned as i */
/* nbin      : bin number */
/* number[i] : number of pairs for bin i */
/* data[i][] : difference of values for bin distance i */
/* v         : sample variogram */
/* Dimension of data is the sum over numbin of number[i] */
{
 int i,j;
 float s;
 for (i = 0; i < nbin; i++)
 {
  s = 0.0;
  for (j = 0; j < number[i]; j++) s += pow(data[i][j],2.0);
  v[i] = 0.5*s/number[i];
 }
 return;
}

void smadogram(int nbin, int *number, float **data, float *v)
/* Isotropic sample madogram calculation */
/* The discrete distance between pair of locations is binned as i */
/* nbin      : bin number */
/* number[i] : number of pairs for bin i */
/* data[i][] : difference of values for bin distance i */
/* v         : sample variogram */
/* Dimension of data is the sum over numbin of number[i] */
{
 int i,j;
 float s;
 for (i = 0; i < nbin; i++)
 {
  s = 0.0;
  for (j = 0; j < number[i]; j++) s += fabs(data[i][j]);
  v[i] = 0.5*s/number[i];
 }
 return;
}
