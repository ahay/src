/* 1D median filtering. */
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

static void extenddata(float* tempt,float* extendt,int nfw,int n1,int n2);
/*extend seismic data*/

int main (int argc, char* argv[]) 
{
	int n1,n2; /*n1 is trace length, n2 is the number of traces*/
	int i;
        int nfw;    /*nfw is the filter-window length*/
	int m;

	float *trace;
	float *tempt;
        float *temp1;
	float *extendt;
    sf_file in, out;

    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);
	/* get the trace length (n1) and the number of traces (n2)*/

    if (!sf_getint("nfw",&nfw)) sf_error("Need integer input");
    /* filter-window length (positive and odd integer)*/
    if (nfw < 1)  sf_error("Need positive integer input"); 
    if (nfw%2 != 0)  nfw = (nfw+1);
    m=(nfw-1)/2;

    trace = sf_floatalloc(n1*n2);
	tempt = sf_floatalloc(n1*n2);
       temp1 = sf_floatalloc(nfw);
	extendt = sf_floatalloc((n1+2*m)*(n2+2*m));

	sf_floatread(trace,n1*n2,in);

	for(int i=0;i<n1*n2;i++)
	{
		tempt[i]=trace[i];
	}

	extenddata(tempt,extendt,nfw,n1,n2);

	   /************1D median filter****************/

	for(i=0;i<n2;i++)
	{
		for(int j=0;j<n1;j++)
		{
			for(int k=0;k<nfw;k++)
			{
				temp1[k]=extendt[(n1+2*m)*i+j+k];
			}
			trace[n1*i+j]=medianfilter(temp1,nfw);
		}
       }

	sf_floatwrite(trace,n1*n2,out);

    exit (0);
}

static void extenddata(float* tempt,float* extendt,int nfw,int n1,int n2)
/*extend seismic data*/
{
    int m=(nfw-1)/2;
    int i;
	for(i=0;i<(n1+2*m)*(n2+2*m);i++)
	{
		extendt[i]=0.0;
	}
	/*extend trace*/
	for(i=0;i<m;i++)
	{
		for(int j=0;j<n1;j++)
		{
			extendt[(n1+2*m)*i+j+m]=0.0;
		}
	}
	for(i=0;i<n2;i++)
	{
		for(int j=0;j<n1;j++)
		{
			extendt[(n1+2*m)*(i+m)+j+m]=tempt[n1*i+j];
		}
	}
	for(i=0;i<m;i++)
	{
		for(int j=0;j<n1;j++)
		{
			extendt[(n1+2*m)*(i+m+n2)+j+m]=0.0;
		}
	}
	/*extend the number of samples*/
	for(i=0;i<(n2+2*m);i++)
	{
		for(int j=0;j<m;j++)
		{
			extendt[(n1+2*m)*i+j]=0.0;
		}
	}
	for(i=0;i<(n2+2*m);i++)
	{
		for(int j=0;j<m;j++)
		{
			extendt[(n1+2*m)*i+j+n1+m]=0.0;
		}
	}
	
}


/* 	$Id: Mmf.c 1131 2007-10-09 23:25:10Z yang $	 */


