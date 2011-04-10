/* 2D Multistage weighted median filtering. */
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

#include "median.h"
#include "weightmf.h"
#include "boundary.h"

int main (int argc, char* argv[]) 
{
    int n1,n2; /*n1 is trace length, n2 is the number of traces*/
    int wn1,wn2;/*n1 is weight data traces length, n2 is the number of weight data traces*/
    int i,j,k,pass,jj;
    int nfw;    /*nfw is the filter-window length*/
    int m;
    float a;   /*temporary variable*/
    bool boundary;  /*boundary parameter*/
    float *trace,*wei;
    float *tempt,*tempw;
    float *extendt,*extendw;
    float *temp1,*temp2,*temp3,*temp4;
    float *tempw1,*tempw2,*tempw3,*tempw4;
    float *z,*Y;
    sf_file in, out, weights;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");
    weights=sf_input("weights");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);
    /* get the trace length (n1) and the number of traces (n2)*/
    
    if (!sf_histint(weights,"n1",&wn1)) sf_error("No n1= in weights");
    wn2 = sf_leftsize(weights,1);
    
    if (n1!=wn1 || n2!=wn2) sf_error("Different size in two input data");
    
    if (!sf_getint("nfw",&nfw)) sf_error("Need integer input");
    /* filter-window length (positive and odd integer)*/

    if (!sf_getbool("boundary",&boundary)) boundary=false;
    /* if y, boundary is data, whereas zero*/

    if (nfw < 1)  sf_error("Need positive integer input"); 
    if (nfw%2 == 0)  nfw = (nfw+1);
    m=(nfw-1)/2;
    
    trace = sf_floatalloc(n1*n2);
    wei = sf_floatalloc(n1*n2);
    tempt = sf_floatalloc(n1*n2);
    extendt = sf_floatalloc((n1+2*m)*(n2+2*m));
    tempw = sf_floatalloc(n1*n2);
    extendw = sf_floatalloc((n1+2*m)*(n2+2*m));
    temp1 = sf_floatalloc(nfw);
    temp2 = sf_floatalloc(nfw);
    temp3 = sf_floatalloc(nfw);
    temp4 = sf_floatalloc(nfw);
    tempw1 = sf_floatalloc(nfw);
    tempw2 = sf_floatalloc(nfw);
    tempw3 = sf_floatalloc(nfw);
    tempw4 = sf_floatalloc(nfw);
    
    z = sf_floatalloc(4);
    Y = sf_floatalloc(3);
    sf_floatread(trace,n1*n2,in);
    sf_floatread(wei,n1*n2,weights);
    
    for(i=0;i<n1*n2;i++)
    {
	tempt[i]=trace[i];
	tempw[i]=fabs(wei[i]);
    }
    
    bound3(tempt,extendt,nfw,nfw,n1,n2,boundary);
    bound3(tempw,extendw,nfw,nfw,n1,n2,boundary);
    
    /************2D multistage weighted median filter****************/
    for(i=m;i<(n2+m);i++){
	for(j=m;j<(n1+m);j++){
	    for(k=0;k<nfw;k++){
		temp1[k]=extendt[(n1+2*m)*i+(j-m+k)];     /*vertical*/
		temp2[k]=extendt[(n1+2*m)*(i-m+k)+j];    /*horizontal*/
		temp3[k]=extendt[(n1+2*m)*(i-m+k)+(j-m+k)];  /*left-up to right-down*/
		temp4[k]=extendt[(n1+2*m)*(i-m+k)+(j+m-k)];   /*left-down to right-up*/
		tempw1[k]=extendw[(n1+2*m)*i+(j-m+k)];     /*vertical*/
		tempw2[k]=extendw[(n1+2*m)*(i-m+k)+j];    /*horizontal*/
		tempw3[k]=extendw[(n1+2*m)*(i-m+k)+(j-m+k)];  /*left-up to right-down*/
		tempw4[k]=extendw[(n1+2*m)*(i-m+k)+(j+m-k)];   /*left-down to right-up*/
	    }
	    z[0]=wmedianfilter(temp1,tempw1,nfw);
	    z[1]=wmedianfilter(temp2,tempw2,nfw);
	    z[2]=wmedianfilter(temp3,tempw3,nfw);
	    z[3]=wmedianfilter(temp4,tempw4,nfw);
	    
	    for(pass=1;pass<4;pass++){
		for(jj=0;jj<4-pass;jj++){
		    if(z[jj]>z[jj+1]){
			a=z[jj];
			z[jj]=z[jj+1];
			z[jj+1]=a;
		    }
		}
	    }
	    Y[0]=z[3];
	    Y[1]=z[0];
	    Y[2]=extendt[(n1+2*m)*i+j];
	    trace[n1*(i-m)+j-m]=medianfilter(Y,3);
	}
    }
    
    sf_floatwrite(trace,n1*n2,out);

    exit (0);
}


/* 	$Id$	 */


