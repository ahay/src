/* Calculate dominant frequency of amplitude spectra dataset.*/
/*
  Copyright (C) 2008 University of Texas at Austin
  
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
    int n1,n2,n3;
    float d1, o1;
    int i,j,k;
    
    float *trace;
    float numer,denom,domin;
    sf_file in;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    if (!sf_histint(in,"n3",&n3)) n3=1;
    /* get the trace length (n1) and the number of traces (n2) and n3*/

    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
   
    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");

   /*set the data space*/
    trace = sf_floatalloc(n1*n2);
    
    for(k = 0; k < n3; k ++) {
	sf_warning("slice %d of %d",(k+1),n3);
	printf("******************************************************\n");
	sf_floatread(trace,n1*n2,in);
	for(j = 0; j < n2; j ++) {
	    numer = 0.;
	    denom = 0.;
	    for(i= 0 ;i < n1; i ++) {
		numer += trace[n1*j+i]*(o1+d1*i);
		denom += trace[n1*j+i];
	    }
	    domin = numer/denom;
	    printf("Dominant frequency at trace  %d    =   %f Hz\n", (j+1),domin);          
	}
    }
    printf("******************************************************\n");

    exit (0);
}

/* 	$Id$	 */
