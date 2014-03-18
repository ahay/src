/* 1D non-local median filtering. */
/*
  Copyright (C) 2014 University of Texas at Austin
  
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
void make_kernal(int f, float *kernal);

int main (int argc, char* argv[]) 
{
    int n1,n2,n3; /*n1 is trace length, n2 is the number of traces, n3 is the number of 3th axis*/
    int nfw;    /*nfw is the filter-window length*/
    int i,j,k,ii;
    int rmin, rmax, t, f, r;
    float wmax,average,sweight,w,d,h;
    bool boundary;
  

    float *trace;
    float *tempt;
    float *w1,*w2,*kernal;
    float *extendt;
    sf_file in, out;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    /* get the trace length (n1) and the number of traces (n2) and n3*/

    if (!sf_getint("t",&t)) t=5;
    /* radio of search window */

    if (!sf_getint("f",&f)) f=2;
    /* radio of similarity window */
	nfw=2*f+1;

    if (!sf_getfloat("h",&h)) h=0.5;
    /* degree of filtering */   
	 
    if (!sf_getbool("boundary",&boundary)) boundary=false;
    /* if y, boundary is data, whereas zero*/
    
    trace = sf_floatalloc(n1*n2);
    memset(trace,0,n1*n2*sizeof(float));

    tempt = sf_floatalloc(n1*n2);
    w1 = sf_floatalloc(nfw);
    w2 = sf_floatalloc(nfw);
    kernal = sf_floatalloc(nfw);
    extendt = sf_floatalloc((n1+2*f)*n2);

    for(ii=0;ii<n3;ii++){
	sf_floatread(trace,n1*n2,in);
	
	for(i=0;i<n1*n2;i++){
	    tempt[i]=trace[i];
	}

	bound1(tempt,extendt,nfw,n1,n2,boundary);
	
	/************1D NLM ****************/

	for(i=0;i<n2;i++)
	{
		for(j=0;j<n1;j++)
		{
			for(k=0;k<nfw;k++){
		   		w1[k]=extendt[(n1+2*f)*i+j+k];
			}

			wmax=0.0;
			average=0.0;
			sweight=0.0;
			rmin=SF_MAX(j-t,0);
			rmax=SF_MIN(j+t,n1-1);
			d=0;
			for(r=rmin;r<=rmax;r++)
			{	
				if(r==i) continue;
				for(k=0;k<nfw;k++){
		   			w2[k]=extendt[(n1+2*f)*i+r+k];
				}	
				make_kernal(f,kernal);		
				for(k=0;k<nfw;k++){
		   			d+=kernal[k]*(w1[k]-w2[k])*(w1[k]-w2[k]);
				}
				w=expf(-d/h);                
                	if(w>wmax) wmax=w;  	
                	sweight += w;
                	average += w*extendt[(n1+2*f)*i+r+f]; 
			}
       		 	average += wmax*extendt[(n1+2*f)*i+j+f];
        			sweight += wmax;
				if(sweight>0) trace[n1*i+j] = average / sweight;
				else trace[n1*i+j] = extendt[(n1+2*f)*i+j+f];
		}
	}

	sf_floatwrite(trace,n1*n2,out);
    }

    exit (0);
}

void make_kernal(int f, float *kernal)
/*< make kernal for NLM >*/
{
    	int d,i,nfw;
    	float value;
	nfw=2*f+1;
     memset(kernal, 0, nfw*sizeof(float));
	for(d=0;d<f;d++)
	{
		value=1.0/(2*d+1)/(2*d+1);
		for(i=0;i<nfw;i++)
			kernal[f+1-i-d] += value;
	}
     for(i=0;i<nfw;i++)
		kernal[i] = kernal[i]/f;
}

