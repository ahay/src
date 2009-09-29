/* 1-D and 2-D modified-trimmed-mean (MTM) filtering. 
Also called range-trimmed-mean filter
1D filter (nfw2=1); 2D filter (otherwise)
median filter (pclip=0.); mean filter (pclip=100.)
*/
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

#include <rsf.h>
#include <stdio.h>
#include <math.h>

#include "mean.h"
#include "boundary.h"

int main (int argc, char* argv[]) 
{
    int n1,n2,n3; /*n1 is trace length, n2 is the number of traces, n3 is the number of 3th axis*/

    int nfw1;    /*nfw is the filter-window length in sample direction*/
    int m1;
    int nfw2;    /*nfw is the filter-window length in trace direction*/
    int m2;
    float pclip;
    char *type;
   
    int i,j,k,ii,kk,jj;
    bool boundary, verb;
    
    float *trace;
    float *tempt;
    float *temp1;
    float *extendt;
    sf_file in, out;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    /* get the trace length (n1) and the number of traces (n2) and n3*/
    
    if (!sf_getint("nfw1",&nfw1)) sf_error("Need integer input");
    /* filter-window length in n1 direction (positive and odd integer)*/
    
    if (!sf_getint("nfw2",&nfw2)) nfw2=1;
    /* filter-window length in n2 direction (default=1, 1D case)*/
    
    if (nfw1 < 1 || nfw2 < 1)  sf_error("Need positive integer input"); 
    if (nfw1%2 == 0)  nfw1 = (nfw1+1);
    if (nfw2%2 == 0)  nfw2 = (nfw2+1);
    m1=(nfw1-1)/2;
    m2=(nfw2-1)/2;
    
    if (!sf_getfloat("pclip",&pclip)) sf_error("Need pclip=");
    /* 0.0 <= pclip <= 100.0: median filter (pclip=0.); mean filter (pclip=100.) */
    if (pclip < 0. || pclip > 100.)  sf_error("Need input is in range [0.0,100.0]"); 

    if (!sf_getbool("boundary",&boundary)) boundary=false;
    /* if y, boundary is data, whereas zero*/

    if (NULL == (type=sf_getstring("type"))) type="rectangular";
    /* [rectangular,cross] 2-D window type, the default is rectangular  */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    trace = sf_floatalloc(n1*n2);
    tempt = sf_floatalloc(n1*n2);
    extendt = sf_floatalloc((n1+2*m1)*(n2+2*m2));

    switch(type[0]) {
	case 'r':
	    temp1 = sf_floatalloc(nfw1*nfw2);
	    for(ii=0;ii<n3;ii++) {
		if (verb) sf_warning("slice %d of %d",ii+1,n3);
		sf_floatread(trace,n1*n2,in);
		for(i=0;i<n1*n2;i++){
		    tempt[i]=trace[i];
		}

		bound3(tempt,extendt,nfw1,nfw2,n1,n2,boundary);

		/************modified-trimmied-mean filter****************/
		for(i=0;i<n2;i++){
		    for(j=0;j<n1;j++) {

			for(k=0;k<nfw2;k++){
			    for(kk=0;kk<nfw1;kk++) {
				temp1[k*nfw1+kk]=extendt[(n1+2*m1)*(i+k)+j+kk];
			    }
			}
			/*square window choosing*/
			trace[n1*i+j]=rangemean(temp1,nfw1*nfw2,pclip);
		    }
		}
		sf_floatwrite(trace,n1*n2,out);
	    }
	    break;
	case 'c':
	    temp1 = sf_floatalloc(nfw1+nfw2-1);
	    for(ii=0;ii<n3;ii++) {
		if (verb) sf_warning("slice %d of %d",ii+1,n3);
		sf_floatread(trace,n1*n2,in);
		for(i=0;i<n1*n2;i++){
		    tempt[i]=trace[i];
		}

		bound3(tempt,extendt,nfw1,nfw2,n1,n2,boundary);

		/************modified-trimmied-mean filter****************/
		for(i=0;i<n2;i++){
		    for(j=0;j<n1;j++) {

			jj=0;
			for(k=0;k<nfw2;k++){
			    for(kk=0;kk<nfw1;kk++) {
				if(k==m2 || kk==m1) {
				    temp1[jj]=extendt[(n1+2*m1)*(i+k)+j+kk];
				    jj++;
				}
			    }
			}
			if(jj!=(nfw1+nfw2-1)) {
			    sf_error("Window calculation error, number(%d)!=nfw1+nwf2-1(%d)",jj,(nfw1+nfw2-1));
			}
			trace[n1*i+j]=rangemean(temp1,nfw1+nfw2-1,pclip);
		    }
		}
		sf_floatwrite(trace,n1*n2,out);
	    }
	    break;
	default:
	    sf_error("Unknown method type=%c",type[0]);
	    break;
    }

    exit (0);
}

/* 	$Id$	 */

