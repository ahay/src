/* 2D non-local median filtering. */
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
void make_kernal(int f1, int f2, float *kernal);

int main (int argc, char* argv[]) 
{
    int n1,n2,n3; /*n1 is trace length, n2 is the number of traces, n3 is the number of 3th axis*/
    int nfw1,nfw2;    /*nfw is the filter-window length*/
    int i,j,k,k1,k2,ii;
    int rmin, rmax,smin,smax, t1,t2,f1,f2, r,s;
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

    if (!sf_getint("t1",&t1)) t1=5;
    /* radio of search window */

    if (!sf_getint("t2",&t2)) t2=5;
    /* radio of search window */

    if (!sf_getint("f1",&f1)) f1=2;
    /* radio of similarity window */

    if (!sf_getint("f2",&f2)) f2=2;
    /* radio of similarity window */
	nfw1=2*f1+1;
	nfw2=2*f2+1;

    if (!sf_getfloat("h",&h)) h=0.5;
    /* degree of filtering */   
	 
    if (!sf_getbool("boundary",&boundary)) boundary=false;
    /* if y, boundary is data, whereas zero*/
    
    trace = sf_floatalloc(n1*n2);
    memset(trace,0,n1*n2*sizeof(float));

    tempt = sf_floatalloc(n1*n2);
    w1 = sf_floatalloc(nfw1*nfw2);
    w2 = sf_floatalloc(nfw1*nfw2);
    kernal = sf_floatalloc(nfw1*nfw2);
    extendt = sf_floatalloc((n1+2*f1)*(n2+2*f2));

    for(ii=0;ii<n3;ii++){
	sf_floatread(trace,n1*n2,in);
	
	for(i=0;i<n1*n2;i++){
	    tempt[i]=trace[i];
	}

	bound3(tempt,extendt,nfw1,nfw2,n1,n2,boundary);

	/************2D NLM ****************/

	for(i=0;i<n2;i++)
	{
		for(j=0;j<n1;j++)
		{
			for(k2=0;k2<nfw2;k2++){
				for(k1=0;k1<nfw1;k1++)
		   			w1[k2*nfw1+k1]=extendt[(n1+2*f1)*(i+k2)+j+k1];
			}

			wmax=0.0;
			average=0.0;
			sweight=0.0;
			rmin=SF_MAX(j-t1,0);
			rmax=SF_MIN(j+t1,n1-1);
			smin=SF_MAX(i-t2,0);
			smax=SF_MIN(i+t2,n2-1);

			d=0;
			for(r=rmin;r<=rmax;r++)
			{	
			for(s=smin;s<=smax;s++)
			{                                               
				if(r==i && s==j) continue;

				for(k2=0;k2<nfw2;k2++){
					for(k1=0;k1<nfw1;k1++)
		   				w2[k2*nfw1+k1]=extendt[(n1+2*f1)*(s+k2)+r+k1];   
				}	

				make_kernal(f1,f2,kernal);	
					
				for(k=0;k<nfw1*nfw2;k++){
		   			d+=kernal[k]*(w1[k]-w2[k])*(w1[k]-w2[k]);
				}

				w=expf(-d/h);                
                	if(w>wmax) wmax=w;  	
                	sweight += w;
                	average += w*extendt[(n1+2*f1)*(s+f2)+r+f1]; 
			}
			}
       		 	average += wmax*extendt[(n1+2*f1)*(i+f2)+j+f1];
        			sweight += wmax;
				if(sweight>0) trace[n1*i+j] = average / sweight;
				else trace[n1*i+j] = extendt[(n1+2*f1)*(i+f2)+j+f1];
		}
	}

	sf_floatwrite(trace,n1*n2,out);
    }

    exit (0);
}

void make_kernal(int f1, int f2, float *kernal)
/*< make kernal for NLM >*/
{
    	int d1,d2,i,j,nfw1,nfw2;
    	float value;
	nfw1=2*f1+1;
	nfw2=2*f2+1;
     memset(kernal, 0, nfw1*nfw2*sizeof(float));

	for(d1=0;d1<f1;d1++)
	{
		value=1.0/(2*d1+1)/(2*d1+1);
		for(i=0;i<nfw2;i++)
			for(j=0;j<nfw1;j++)
				kernal[nfw1*(f2+1-i-d1)+f1+1-j-d1] += value;
	}
     for(i=0;i<nfw1*nfw2;i++)
		kernal[i] = kernal[i]/f1;
}

//void make_kernal(int f1, int f2, float *kernal)
/*< make kernal for NLM >*/
/*{
    	int d1,d2,i,j,nfw1,nfw2;
    	float value;
	nfw1=2*f1+1;
	nfw2=2*f2+1;
     memset(kernal, 0, nfw1*nfw2*sizeof(float));
	for(d2=0;d2<f2;d2++)
	for(d1=0;d1<f1;d1++)
	{
		value=1.0/(2*d1+1)/(2*d2+1);
		for(i=0;i<nfw2;i++)
			for(j=0;j<nfw1;j++)
				kernal[nfw1*(f2+1-i-d2)+f1+1-j-d1] += value;
	}
     for(i=0;i<nfw1*nfw2;i++)
		kernal[i] = kernal[i]/f1/f2;
}*/

