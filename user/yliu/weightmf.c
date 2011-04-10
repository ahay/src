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
    free(weights);
    free(data);
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
    free(index);
    free(data);
    return wm;
}

void cweight(float *data,float *weight,float nw,float rect)
/*< get weights from data using local correlation >*/
{   
    int i, j, m1, m2;
    float *one, *two, done, dtwo;

    m1 = (nw-1)/2;
    m2 = (rect-1)/2;

    one = sf_floatalloc(rect);
    two = sf_floatalloc(rect);

    done = 0.;
    for (i=0; i < rect; i++) {
	one[i] = data[m1-m2+i];
	done += one[i]*one[i];
    }

    for (i=0; i < nw; i++) {
	for (j=0; j < rect; j++) {
	    if ((i-m2+j) >=0 && (i-m2+j) < nw) {
		two[j] = data[i-m2+j];
	    } else {
		two[j] = 0.;
	    }
	}
	dtwo = 0.;
	weight[i] = 0.;
	for (j=0; j < rect; j++) {
	    dtwo += two[j]*two[j];
	    weight[i] += one[j]*two[j];
	}
	weight[i] = fabsf(weight[i])/(sqrtf(done*dtwo)+FLT_EPSILON);
    }
    free(one);
    free(two);
}

void vweight(float *data,float *weight,float nw)
/*< get weights from data using variance >*/
{   
    int i;
    float mean, norm;

    mean = 0.;
    for (i=0; i < nw; i++) {
	mean += data[i];
    }
    mean = mean/nw;
    norm = 0.;
    for (i=0; i < nw; i++) {
	weight[i] = 1./((data[i]-mean)*(data[i]-mean)+FLT_EPSILON);
	norm +=weight[i];
    }
    for (i=0; i < nw; i++) {
	weight[i] = nw*weight[i]/(norm+FLT_EPSILON);
    }  
}

void lvweight(float *data,float *weight,float nw,float rect)
/*< get weights from data using local variance >*/
{   
    int i, j, m1, m2;
    float *one, *two, mean;

    m1 = (nw-1)/2;
    m2 = (rect-1)/2;
    one = sf_floatalloc(rect);   
    two = sf_floatalloc(rect);

    mean = 0.;
    for (i=0; i < rect; i++) {
	one[i] = data[m1-m2+i];
	mean += one[i];
    }
    mean = mean/rect;
    for (i=0; i < nw; i++) {
	for (j=0; j < rect; j++) {
	    if ((i-m2+j) >=0 && (i-m2+j) < nw) {
		two[j] = data[i-m2+j];
	    } else {
		two[j] = 0.;
	    }
	}
	weight[i] = 0.;
	for (j=0; j < rect; j++) {
	    weight[i] += (two[j]-mean)*(two[j]-mean);
	}
	weight[i] = rect*1./(weight[i]+FLT_EPSILON);
    }
    free(one);
    free(two);
}

/* 	$Id$	 */


