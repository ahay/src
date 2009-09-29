/* 1D LUM filter*/
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

#include "lum.h"
#include "boundary.h"

int main (int argc, char* argv[]) 
{
    int n1,n2,n3; /*n1 is trace length, n2 is the number of traces, n3 is the number of 3th axis*/
    int nfw;    /*nfw is the filter-window length*/
    int smnclip;  /*smnclip is smoother tuning parameter(<=(nfw+1)/2)*/
    int shnclip;  /*shnclip is sharpener tuning parameter(<=(nfw+1)/2)*/
    
    int m;
    int i,j,k,ii;
    bool boundary;
    
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
    
    if (!sf_getint("nfw",&nfw)) sf_error("Need integer input");
    /* filter-window length (positive and odd integer)*/
    
    if (nfw < 1)  sf_error("Need positive integer input"); 
    if (nfw%2 == 0)  nfw = (nfw+1);
    m=(nfw-1)/2;
    
    if (!sf_getint("smnclip",&smnclip)) smnclip=(nfw+1)/2;
    /* smoother tuning parameter (1 <= smnclip <= (nfw+1)/2, the default is (nfw+1)/2)*/

    if (!sf_getint("shnclip",&shnclip)) shnclip=(nfw+1)/2;
    /* sharpener tuning parameter (1 <= shnclip <= (nfw+1)/2, the default is (nfw+1)/2)*/

    if (!sf_getbool("boundary",&boundary)) boundary=false;
    /* if y, boundary is data, whereas zero*/

    if(smnclip<1 || smnclip>((nfw+1)/2)) sf_error("Smnclip need a value between 1 and (nfw+1)/2");
    if(shnclip<1 || shnclip>((nfw+1)/2)) sf_error("Shnclip need a value between 1 and (nfw+1)/2");
    if(smnclip>shnclip) sf_error("Smnclip cannot be greater than shnclip");
    
    trace = sf_floatalloc(n1*n2);
    tempt = sf_floatalloc(n1*n2);
    temp1 = sf_floatalloc(nfw);
    extendt = sf_floatalloc((n1+2*m)*n2);
    
    for(ii=0;ii<n3;ii++) {
	sf_floatread(trace,n1*n2,in);
	
	for(i=0;i<n1*n2;i++){
	    tempt[i]=trace[i];
	}
	
	bound1(tempt,extendt,nfw,n1,n2,boundary);
	
	/************1D LUM smoother filter****************/
	
	for(i=0;i<n2;i++){
	    for(j=0;j<n1;j++) {
		for(k=0;k<nfw;k++){
		    temp1[k]=extendt[(n1+2*m)*i+j+k];
		}
		trace[n1*i+j]=lum(temp1,nfw,smnclip,shnclip);
	    }
	}
	
	sf_floatwrite(trace,n1*n2,out);
    }

    exit (0);
}

/* 	$Id$	 */
