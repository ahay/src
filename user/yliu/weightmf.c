/* 1-D weighted median filtering. */
/*
  Copyright (C) 2007 University of Texas at Austin
  
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

#include <stdio.h>
#include <math.h>
#include <rsf.h>

#include "weightmf.h"

float wmedianfilter(float *odata, float *oweights, int nfw)
/*< get a weighted median value >*/
{   
    int i,j,k;
    float *data, *weights;
    float temp,average,wm,sum;
    
    data = sf_floatalloc(nfw);
    weights = sf_floatalloc(nfw);
    
    for(i=0;i<nfw;i++){
	data[i]=odata[i];
	weights[i]=oweights[i];
    }
    
    for(i=1;i<nfw;i++){
	for(j=0;j<nfw-i;j++){
	    if(data[j]>data[j+1]){
		temp=data[j];
		data[j]=data[j+1];
		data[j+1]=temp;
		temp=weights[j];
		weights[j]=weights[j+1];
		weights[j+1]=temp;
	    }
	}
    }
    sum=0.;
    for(i=0;i<nfw;i++){
	sum+=weights[i];
    }
    average=sum/2.;
    k=0;
    sum=weights[0];
    while (sum<average){
	k++;
	sum+=weights[k];
    }
    wm=data[k];
    return wm;
    
}

float wmedian(float *temp,float *weight,int nfw)
/*< get a weighted median value >*/
{   
    int i, j, pass, *index, b;
    float wm, a, *data;
    
    data = sf_floatalloc(nfw);
    index = sf_intalloc(nfw);
    for (j=0; j < nfw; j++) {
	data[j] = temp[j]*weight[j];
	index[j] = j;
    }
    for(pass=1; pass < nfw; pass++) {
	for(i=0; i < nfw-pass; i++) {
	    if(data[i] > data[i+1]) {
		a = data[i];
		b = index[i];
		data[i] = data[i+1];
		index[i] = index[i+1];
		data[i+1] = a;
		index[i+1] = b;
	    }
	}
    }
    wm = temp[index[nfw/2]];
    return wm;
}

/* 	$Id$	 */


