/* Two-step space varying median filtering. 
   In default case, sftsmf is equal to sftvmf.
*/
/*
  Copyright (C) 2013 University of Texas at Austin
  
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

#include "medianutil.h"

int main (int argc, char* argv[]) 
{
    int n1,n2,n3,N; /*n1 is trace length, n2 is the number of traces, n3 is the number of 3th axis, N is for computing average energy level*/
    int i,j,k,kk,ii,f1,f2; /*f1 is the processing window starting point and f2 is the ending point*/
    int nfw;    /*nfw is the reference filter-window length*/
    int tempnfw;  /*temporary variable*/
    int m;
    float medianv,ael; /*temporary median variable*/
    bool boundary,verb;
    int l1,l2,l3,l4; /*space-varying window coefficients*/
    
    float *trace;
    float *tempt; /*temporary array*/
    float *length; /*output variable filter length map*/
    float *result; /*output array*/
    float *extendt;
    float *medianarray;   /*1D median filtered array*/
    float *temp1,*temp2,*temp3; /*temporary array*/
    sf_file in, out, lengthout;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");
    lengthout=sf_output("L");		

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    /* get the trace length (n1) and the number of traces (n2) and n3*/
    
    if (!sf_getbool("boundary",&boundary)) boundary=false;
    /* if y, boundary is data, whereas zero*/

    if (!sf_getint("ns", &f1)) f1=0;
    /* processing window starting point, corresponding to the temporal axis */

    if (!sf_getint("ne", &f2)) f2=n2-1;
    /* processing window ending point, corresponding to the temporal axis, n2 means transposed first-axis dimension. */

    if (!sf_getint("N", &N))   N=(f2-f1+1)*n1;
    /* average energy level (AEL) computing number */

    if (!sf_getfloat("ael",&ael)) ael=0.0;
    /*	get the average energy level (AEL) empirically defined */

    if (!sf_getbool("verb",&verb)) verb=false;
    /*	if print the computed average energy level (AEL) */
	
    if (!sf_getint("nfw",&nfw)) sf_error("Need integer input");
    /* reference filter-window length (>l4, positive and odd integer)*/
    
    if (!sf_getint("l1",&l1)) l1=2;
    /* space-varying window parameter "l1" (default=2)*/
    
    if (!sf_getint("l2",&l2)) l2=0;
    /* space-varying window parameter "l2" (default=0)*/
    
    if (!sf_getint("l3",&l3)) l3=2;
    /* space-varying window parameter "l3" (default=2)*/
    
    if (!sf_getint("l4",&l4)) l4=4;
    /* space-varying window parameter "l4" (default=4)*/

    if(NULL!=sf_getstring("L")) {
	length=sf_floatalloc(n1*n2);
	for(i=0;i<n2;i++)
	    for(j=0;j<n1;j++)
		length[i*n1+j]=0;
    } else {
	length=NULL;
    }

    if (l1<l2 || l4<l3) sf_error("Need l1>=l2 && l4>=l3"); 
    if ((l1%2)!=0) l1 = l1+1;
    if ((l2%2)!=0) l2 = l2+1;
    if ((l3%2)!=0) l3 = l3+1;
    if ((l4%2)!=0) l4 = l4+1;
    
    if (nfw <=l4)  sf_error("Need nfw > l4"); 
    if (nfw%2 == 0)  nfw = (nfw+1);
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
	    result[i]=trace[i];
	}
	
	bound1(tempt,extendt,nfw,n1,n2,boundary);
	
	if(ael==0.0)
	{/************1D reference median filtering****************/
	
	    for(i=f1;i<f2+1;i++){
		for(j=0;j<n1;j++){
		    for(k=0;k<nfw;k++){
			temp1[k]=extendt[(n1+2*m)*i+j+k];
		    }
		    medianarray[n1*i+j]=medianfilter(temp1,nfw);
		}
	    }
	    medianv=0.0;
	    for(i=f1;i<f2+1;i++)
		for(j=0;j<n1;j++)
		    medianv=medianv+fabs(medianarray[n1*i+j]);
	    
	    medianv=medianv/(1.0*N);}
	else
	{medianv=ael;}
	
	if(verb) sf_warning("The average energy level (AEL) is %g",medianv);
	
	/************1D space-varying median filter****************/
	for(i=f1;i<f2+1;i++){
	    for(kk=0;kk<n1;kk++){
		temp2[kk]=trace[n1*i+kk];
	    }
	    for(j=0;j<n1;j++){
		if(fabs(medianarray[n1*i+j])<medianv){
		    if(fabs(medianarray[n1*i+j])<medianv/2.0){
			tempnfw=nfw+l1;
		    }
		    else{
			tempnfw=nfw+l2;
		    }
		}
		else{
		    if(fabs(medianarray[n1*i+j])>=(medianv*2.0)){
			tempnfw=nfw-l4;
		    }
		    else{
			tempnfw=nfw-l3;
		    }
		}
		if(NULL!=sf_getstring("L")) length[n1*i+j]	= tempnfw;		
		temp3 = sf_floatalloc(tempnfw);
		bound2(temp2,temp3,n1,tempnfw,j,boundary);
		result[n1*i+j]=medianfilter(temp3,tempnfw);
		tempnfw=nfw;
	    }
	}
	sf_floatwrite(result,n1*n2,out);
	if(NULL!=sf_getstring("L"))  sf_floatwrite(length,n1*n2,lengthout);	
    }

    exit (0);
}



