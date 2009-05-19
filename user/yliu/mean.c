/* Mean filtering. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

#include "mean.h"

float mean(float *a, int n)
/*< get a mean value >*/
{
    int i;
    float t=0;

    for (i=0; i < n; i++) 
	t+=a[i];
    t = t/n;

    return t;
}

float alphamean(float *temp, int nfw, float alpha)
/*< get a alpha-trimmed-mean value >*/
{
    int i,j,pass,lowc,higc;
    float mean,a;
    for(pass=1;pass<nfw;pass++){
	for(i=0;i<nfw-pass;i++){
	    if(temp[i]>temp[i+1]){
		a=temp[i];
		temp[i]=temp[i+1];
		temp[i+1]=a;
	    }
	}
    }  
    lowc=alpha*nfw;
    higc=(1-alpha)*nfw;
    mean=0.;
    j=lowc;
    do {
	mean+=temp[j];
	j++;
    }while(j<higc);
    mean/=(1-2*alpha)*(nfw-1)+1;
    return mean;
}

float rangemean(float *temp, int nfw, float pclip)
/*< get a modified-trimmed-mean value >*/
{
    int i,pass,n;
    float median,mean,a,dis1,dis2,*data,*weight;
    data = sf_floatalloc(nfw);
    weight = sf_floatalloc(nfw);
    for(pass=0;pass<nfw;pass++){
	data[pass]=temp[pass];
    }
    for(pass=1;pass<nfw;pass++){
	for(i=0;i<nfw-pass;i++){
	    if(temp[i]>temp[i+1]){
		a=temp[i];
		temp[i]=temp[i+1];
		temp[i+1]=a;
	    }
	}
    }  
    median=temp[nfw/2];
    dis1=fabs(temp[0]-median);
    dis2=fabs(temp[nfw-1]-median);
    mean=0.;
    n=0;
    if (dis1 >= dis2){
	pclip*=dis1*0.01;
    } else {
	pclip*=dis2*0.01;
    }
    for(pass=0;pass<nfw;pass++){
	if (data[pass]>=(median-pclip) && data[pass]<=(median+pclip)){
	    weight[pass]=1.;
	    n++;
	} else{
	    weight[pass]=0.;
	}
	mean+=data[pass]*weight[pass];
    }
    if (n==0) n=1;
    mean/=n;
    return mean;
}

/* 	$Id$	 */


