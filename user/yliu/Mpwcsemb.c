/* Semblance from plane-wave construction datacube. */
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

#include "semblance.h"
#include "boundary.h"

int main (int argc, char* argv[]) 
{
    int n1,n2,n3; /*n1 is trace length, n2 is the number of prediction, n3 is the number of trace */

    int nfw;    /*nfw is the calculation window in sample direction*/
    int m;
   
    int i,k,ii,kk;
    bool boundary, verb;
    
    float *trace, *tempt, *temp1, *extendt, *result;
    sf_file in, out;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    sf_unshiftdim(in, out, 2);
    sf_putint(out,"n3",1);

    /* get the trace length (n1) and the number of traces (n2) and n3*/
    
    if (!sf_getint("nfw",&nfw)) sf_error("Need integer input");
    /* calculation window in n1 direction (positive integer)*/
    
    if (nfw < 1)  sf_error("Need positive integer input"); 
    if (nfw%2 == 0)  nfw = (nfw+1);
    m = (nfw-1)/2;
    
    if (!sf_getbool("boundary",&boundary)) boundary=false;
    /* if y, boundary is data, whereas zero*/

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    trace = sf_floatalloc(n1*n2);
    tempt = sf_floatalloc(n1*n2);
    extendt = sf_floatalloc((n1+2*m)*n2);

    temp1 = sf_floatalloc(nfw*n2);
    result = sf_floatalloc(n1);
    for(ii=0;ii<n3;ii++) {
	if (verb) sf_warning("slice %d of %d",ii+1,n3);
	sf_floatread(trace,n1*n2,in);
	for(i=0;i<n1*n2;i++){
	    tempt[i]=trace[i];
	}
	
	bound3(tempt,extendt,nfw,1,n1,n2,boundary);
	
	for(i=0;i<n1;i++) {
	    
	    for(k=0;k<n2;k++){
		for(kk=0;kk<nfw;kk++) {
		    temp1[k*nfw+kk]=extendt[(n1+2*m)*k+i+kk];
		}
	    }
	    result[i]=semblance(temp1,nfw,n2);
	}
	sf_floatwrite(result,n1,out);
    }

    exit (0);
}

/* 	$Id$	 */

