/* 1D Time-varying median filtering. */
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

#include <rsf.h>
#include <stdio.h>
#include <math.h>

#include "boundary.h"

int main (int argc, char* argv[]) 
{
    int n1,n2,n3; /*n1 is trace length, n2 is the number of traces, n3 is the number of 3th axis*/
    int i,j,k,kk,ii;
    int nfw;    /*nfw is the reference filter-window length*/
    int tempnfw;  /*temporary variable*/
    int m;
    float medianv; /*temporary median variable*/
    bool boundary;
    int alpha,beta,gamma,delta; /*time-varying window coefficients*/
    
    float *trace;
    float *tempt; /*temporary array*/
    float *result; /*output array*/
    float *extendt;
    float *medianarray;   /*1D median filtered array*/
    float *temp1,*temp2,*temp3; /*temporary array*/
    sf_file in, out;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    /* get the trace length (n1) and the number of traces (n2) and n3*/
    
    if (!sf_getbool("boundary",&boundary)) boundary=false;
    /* if y, boundary is data, whereas zero*/
    
    if (!sf_getint("nfw",&nfw)) sf_error("Need integer input");
    /* reference filter-window length (>delta, positive and odd integer)*/
    
    if (!sf_getint("alpha",&alpha)) alpha=2;
    /* time-varying window parameter "alpha" (default=2)*/
    
    if (!sf_getint("beta",&beta)) beta=0;
    /* time-varying window parameter "beta" (default=0)*/
    
    if (!sf_getint("gamma",&gamma)) gamma=2;
    /* time-varying window parameter "gamma" (default=2)*/
    
    if (!sf_getint("delta",&delta)) delta=4;
    /* time-varying window parameter "delta" (default=4)*/
    
    if (alpha<beta || delta<gamma) sf_error("Need alpha>=beta && delta>=gamma"); 
    if ((alpha%2)!=0) alpha = alpha+1;
    if ((beta%2)!=0) beta = beta+1;
    if ((gamma%2)!=0) gamma = gamma+1;
    if ((delta%2)!=0) delta = delta+1;
    
    if (nfw <=delta)  sf_error("Need nfw > delta"); 
    if (nfw%2 == 0)  nfw++;
    m=(nfw-1)/2;
    tempnfw=nfw;
    
    trace = sf_floatalloc(n1*n2);
    tempt = sf_floatalloc(n1*n2);
    result = sf_floatalloc(n1*n2);
    extendt = sf_floatalloc((n1+2*m)*n2);
    medianarray =sf_floatalloc(n1*n2);
    temp1 = sf_floatalloc(nfw);
    temp2 = sf_floatalloc(n1);
    /*set the data space*/
    
    for(ii=0;ii<n3;ii++){
	sf_floatread(trace,n1*n2,in);
	for(i=0;i<n1*n2;i++){
	    tempt[i]=trace[i];
	}
	
	bound1(tempt,extendt,nfw,n1,n2,boundary);
	
	/************1D reference median filtering****************/
	
	for(i=0;i<n2;i++){
	    for(j=0;j<n1;j++){
		for(k=0;k<nfw;k++){
		    temp1[k]=extendt[(n1+2*m)*i+j+k];
		}
		medianarray[n1*i+j]=sf_quantile(m,nfw,temp1); 
	    }
	}
	medianv=0.0;
	for(i=0;i<n1*n2;i++){
	    medianv=medianv+fabs(medianarray[i]);
	}
	medianv=medianv/(1.0*n1*n2);
	
	/************1D time-varying median filter****************/
	for(i=0;i<n2;i++){
	    for(kk=0;kk<n1;kk++){
		temp2[kk]=trace[n1*i+kk];
	    }
	    for(j=0;j<n1;j++){
		if(fabs(medianarray[n1*i+j])<medianv){
		    if(fabs(medianarray[n1*i+j])<medianv/2.0){
			tempnfw=nfw+alpha;
		    }
		    else{
			tempnfw=nfw+beta;
		    }
		}
		else{
		    if(fabs(medianarray[n1*i+j])>=(medianv*2.0)){
			tempnfw=nfw-delta;
		    }
		    else{
			tempnfw=nfw-gamma;
		    }
		}
		temp3 = sf_floatalloc(tempnfw);
		bound2(temp2,temp3,n1,tempnfw,j,boundary);
		result[n1*i+j]=sf_quantile((tempnfw-1)/2,tempnfw,temp3); 
		tempnfw=nfw;
	    }
	}
	sf_floatwrite(result,n1*n2,out);
    }

    exit (0);
}


/* 	$Id: Mtvmf.c 3368 2008-03-09 20:42:10Z yang	 */


