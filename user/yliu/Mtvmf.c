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

#include "median.h"

void extenddata1(float* tempt,float* extendt,int nfw,int n1,int n2);/*extend seismic data*/
void extenddata2(float* temp2,float* temp3,int n1,int tempnfw,int j);/*extend temporary seismic data*/
float medianfilter(float* temp,int nfw); /*get the median value from a queue*/

int main (int argc, char* argv[]) 
{
	int n1,n2; /*n1 is trace length, n2 is the number of traces*/
	int i,j,k,kk;
	int nfw;    /*nfw is the reference filter-window length*/
	int tempnfw;  /*temporary variable*/
	int m;
	float medianv; /*temporary median variable*/


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
	n2 = sf_leftsize(in,1);
	/* get the trace length (n1) and the number of traces (n2)*/

	if (!sf_getint("nfw",&nfw)) sf_error("Need integer input");
	/* reference filter-window length (>7, positive and odd integer)*/
	if (nfw < 7)  sf_error("Need positive integer input and greater than 7"); 
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

	sf_floatread(trace,n1*n2,in);

	for(i=0;i<n1*n2;i++)
	{
		tempt[i]=trace[i];
	}

	extenddata1(tempt,extendt,nfw,n1,n2);

	/************1D reference median filtering****************/

	for(i=0;i<n2;i++)
	{
		for(j=0;j<n1;j++)
		{
			for(k=0;k<nfw;k++)
			{
				temp1[k]=extendt[(n1+2*m)*i+j+k];
			}
			medianarray[n1*i+j]=medianfilter(temp1,nfw);
		}
       }
	medianv=0.0;
	for(i=0;i<n1*n2;i++)
	{
		medianv=medianv+fabs(medianarray[i]);
	}
	medianv=medianv/(1.0*n1*n2);

	/************1D time-varying median filter****************/
	for(i=0;i<n2;i++)
	{
		for(kk=0;kk<n1;kk++)
		{
			temp2[kk]=trace[n1*i+kk];
		}
		for(j=0;j<n1;j++)
		{
			if(fabs(medianarray[n1*i+j])<medianv)
			{
				if(fabs(medianarray[n1*i+j])<medianv/2.0)
				{
					tempnfw=nfw+2;
				}
				else
				{
					tempnfw=nfw;
				}
			}
			else
			{
				if(fabs(medianarray[n1*i+j])>=(medianv*2.0))
				{
					tempnfw=nfw-4;
				}
				else
				{
					tempnfw=nfw-2;
				}
			}
			temp3 = sf_floatalloc(tempnfw);
			extenddata2(temp2,temp3,n1,tempnfw,j);
			result[n1*i+j]=medianfilter(temp3,tempnfw);
			tempnfw=nfw;
		
		}
       }

	sf_floatwrite(result,n1*n2,out);

    exit (0);
}

void extenddata1(float* tempt,float* extendt,int nfw,int n1,int n2)/*extend seismic data*/
{
	int m=(nfw-1)/2;
	int i,j;

	for(i=0;i<(n1+2*m)*(n2);i++)
	{
		extendt[i]=0.0;
	}
	/*extend the number of samples*/
	for(i=0;i<n2;i++)
	{
		for(j=0;j<m;j++)
		{
			extendt[(n1+2*m)*i+j]=0.0;
		}
	}
	for(i=0;i<n2;i++)
	{
		for(j=0;j<n1;j++)
		{
			extendt[(n1+2*m)*i+j+m]=tempt[n1*i+j];
		}
	}
	for(i=0;i<n2;i++)
	{
		for(j=0;j<m;j++)
		{
			extendt[(n1+2*m)*i+j+n1+m]=0.0;
		}
	}
}

void extenddata2(float* temp2,float* temp3,int n1,int tempnfw,int j)/*extend temporary seismic data*/
{
	int k;
	/*extend trace*/
	if((j-tempnfw/2)>=0&&(j+tempnfw/2)<n1)
	{
		for(k=0;k<tempnfw;k++)
		{
			temp3[k]=temp2[(j-tempnfw/2+k)];
		}
	}
	else if((j-tempnfw/2)<0&&(j+tempnfw/2)<n1)
	{
		for(k=0;k<(abs(j-tempnfw/2));k++)
		{
			temp3[k]=0.0;
		}
		for(k=(abs(j-tempnfw/2));k<tempnfw;k++)
		{
			temp3[k]=temp2[k-abs(j-tempnfw/2)];
		}
	}
	else if((j-tempnfw/2)>=0&&(j+tempnfw/2)>=n1)
	{
		for(k=0;k<(tempnfw-abs(j+tempnfw/2-n1+1));k++)
		{
			temp3[k]=temp2[j-tempnfw/2+k];
		}
		for(k=(tempnfw-abs(j+tempnfw/2-n1+1));k<tempnfw;k++)
		{
			temp3[k]=0.0;
		}
	}
	else
	{
		for(k=0;k<(abs(j-tempnfw/2));k++)
		{
			temp3[k]=0.0;
		}
		for(k=(abs(j-tempnfw/2));k<(tempnfw-abs(j+tempnfw/2-n1+1));k++)
		{
			temp3[k]=temp2[k-abs(j-tempnfw/2)];
		}
		for(k=(tempnfw-abs(j+tempnfw/2-n1+1));k<tempnfw;k++)
		{
			temp3[k]=0.0;
		}
	}	
}

/* 	$Id: Mtvmf.c 3303 2008-02-06 15:17:10Z yang $	 */


