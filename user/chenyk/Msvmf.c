/* Space varying median filtering. 
   Using local similarity as a reference.
*/
/*
  Copyright (C) 2014 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as publishIn default case, sftsmf is equal to sftvmf.ed by
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
    int n1,n2,n3; /*n1 is trace length, n2 is the number of traces, n3 is the number of 3th axis, N is for computing average energy level*/
    int i,j,kk,ii,f1,f2; /*f1 is the processing window starting point and f2 is the ending point*/
    int nfw;    /*nfw is the reference filter-window length*/
    int tempnfw;  /*temporary variable*/
    bool boundary;
    int l1,l2,l3,l4; /*space-varying window coefficients*/
    float lambda1,lambda2,lambda3,lambda4; /*space-varying window coefficients*/   
     
    float *trace;
    float *simi; /*similarity array*/
    float *simi1; /*similarity array*/
    float s1,s2,s3,s4;
    float *length; /*output variable filter length map*/
    float *result; /*output array*/

    float *temp2,*temp3; /*temporary array*/
    sf_file in, out, lengthout,similarity;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");
    similarity = sf_input("similarity");
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
	
    if (!sf_getint("nfw",&nfw)) sf_error("Need integer input");
    /* reference filter-window length (>l4, positive and odd integer)*/
    
    if (!sf_getint("l1",&l1)) l1=4;
    /* space-varying window parameter "l1" (default=4)*/
    
    if (!sf_getint("l2",&l2)) l2=2;
    /* space-varying window parameter "l2" (default=2)*/
    
    if (!sf_getint("l3",&l3)) l3=2;
    /* space-varying window parameter "l3" (default=2)*/
    
    if (!sf_getint("l4",&l4)) l4=4;
    /* space-varying window parameter "l4" (default=4)*/


    if (!sf_getfloat("lambda1",&lambda1)) lambda1=0.15;
    /* space-varying window parameter "lambda1" (default=0.15)*/
    
    if (!sf_getfloat("lambda2",&lambda2)) lambda2=0.25;
    /* space-varying window parameter "lambda2" (default=0.25)*/
    
    if (!sf_getfloat("lambda3",&lambda3)) lambda3=0.75;
    /* space-varying window parameter "lambda3" (default=0.75)*/
    
    if (!sf_getfloat("lambda4",&lambda4)) lambda4=0.85;
    /* space-varying window parameter "lambda4" (default=0.85)*/
    
    
    if (l1<l2 || l4<l3) sf_error("Need l1>=l2 && l4>=l3"); 
    if ((l1%2)!=0) l1 = l1+1;
    if ((l2%2)!=0) l2 = l2+1;
    if ((l3%2)!=0) l3 = l3+1;
    if ((l4%2)!=0) l4 = l4+1;
    
    if (nfw <=l4)  sf_error("Need nfw > l4"); 
    if (nfw%2 == 0)  nfw = (nfw+1);
    tempnfw=nfw;
    
    trace = sf_floatalloc(n1*n2);
    simi = sf_floatalloc(n1*n2);
    simi1 = sf_floatalloc(n1*n2);
    result = sf_floatalloc(n1*n2);

    temp2 = sf_floatalloc(n1);
    /*set the data space*/
    
    if(NULL!=sf_getstring("L")) {
	length=sf_floatalloc(n1*n2);
    } else {
	length=NULL;
    }


    for(ii=0;ii<n3;ii++){

	if(NULL!=length)
		{
		for(i=0;i<n2;i++)
			for(j=0;j<n1;j++)
				length[i*n1+j]=0;
		}

	sf_floatread(trace,n1*n2,in);
	sf_floatread(simi,n1*n2,similarity);

	for(i=0;i<n1*n2;i++){
	    result[i]=0;
		simi1[i]=simi[i];
	}
	
	s1=   sf_quantile((int)(lambda1*n1*n2),n1*n2,simi1); 	
	s2=   sf_quantile((int)(lambda2*n1*n2),n1*n2,simi1); 	
	s3=   sf_quantile((int)(lambda3*n1*n2),n1*n2,simi1); 	
	s4=   sf_quantile((int)(lambda4*n1*n2),n1*n2,simi1); 	

	/************1D space-varying median filter****************/
	for(i=f1;i<f2+1;i++){
	    for(kk=0;kk<n1;kk++){
		temp2[kk]=trace[n1*i+kk];
	    }
	    for(j=0;j<n1;j++){
		if(simi[n1*i+j]<s2){
		    if(simi[n1*i+j]<=s1){
			tempnfw=nfw+l1;
			/* sf_warning("simi[%d*%d+%d]=%g,L=%d",n1,i,j,simi[n1*i+j],tempnfw); */
		    }
		    else{
			tempnfw=nfw+l2;
		    }
		}
		else {
			if(simi[n1*i+j]>s3)
			{
		    		if(simi[n1*i+j]>s4){
					tempnfw=nfw-l4;
		    		}
		    		else{
					tempnfw=nfw-l3;
		    		}
			}
		}
		if(NULL!=length) length[n1*i+j]	= tempnfw;		
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



