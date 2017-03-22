/* Instantaneous signal-to-noise ratio (SNR) estimation. */
/*
  Copyright (C) 2017 Jilin University
  
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

int main (int argc, char* argv[]) 
{
    int n1,n2,n3,n12; /*n1 is trace length, n2 is the number of traces*/
    int i,j,k,l,m,h;
    float *s, *n, *output, max, min;

    sf_file in, out, en;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out= sf_output("out");
    en = sf_input("en");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(en,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(en,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    n12 = n1*n2;
    
    s = sf_floatalloc(n12);
    output = sf_floatalloc(n12);
    n = sf_floatalloc(n12);
    
    sf_floatread(s,n12,in);
    sf_floatread(n,n12,en);
    
    for(i=0;i<n1;i++) {
	for(j=0;j<n2;j++) {
	    if(fabs(s[n2*i+j]) <= FLT_EPSILON){
		output[n2*i+j]=0;}
	    else{
		if(fabs(n[n2*i+j]) <= FLT_EPSILON){
		    output[n2*i+j]=0;}
		
		else{
		    output[n2*i+j]=10*log10((s[n2*i+j]*s[n2*i+j])/
					    (n[n2*i+j]*n[n2*i+j]));
		}
	    }
	}
    }
    for(k=0;k<n1;k++) {
	for(l=0;l<n2;l++) {
	    max=output[0];
	    if(max<output[k*n2+l]){
		max=output[k*n2+l];
	    }
	    min=output[0];
	    if(min>output[k*n2+l]){
		min=output[k*n2+l];
	    }
	    
	    for(m=0;m<n1;m++) {
		for(h=0;h<n2;h++) {
		    if(fabs(s[n2*m+h]) <= FLT_EPSILON){
			output[n2*m+h]=min;
		    }
		    else{
			if(fabs(n[n2*m+h]) <= FLT_EPSILON){
			    output[n2*m+h]=max;
			}
		    }
		}
	    }
	}
    }  
    sf_floatwrite(output,n12,out);
    
    exit (0);
    
}
/* 	$Id$	 */

