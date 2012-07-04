/* 1D denoising using edge-preserving local polynomial fitting (ELPF). */
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

#include <math.h>
#include "elpf.h"

int main (int argc, char* argv[]) 
{
    int n1, n2, nfw, i2, rect, niter;
    bool boundary, verb;
    
    float *input, *output;
    sf_file in, out;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);
    /* get the trace length (n1) and the number of traces (n2) and n3*/
    
    if (!sf_getint("nfw",&nfw)) sf_error("Need integer input");
    /* filter-window length (positive and odd integer) */

    if (!sf_getint("rect",&rect)) sf_error("Need integer input");
    /* local smoothing radius */

    if (!sf_getbool("boundary",&boundary)) boundary = false;
    /* if y, boundary is data, whereas zero*/

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    
    input = sf_floatalloc(n1);
    output = sf_floatalloc(n1);

    elpf_init(n1,nfw,rect,verb);
   
    for(i2=0; i2 < n2; i2++) {
	sf_floatread(input,n1,in);
	
	elpf (input,output,n1,nfw,niter,boundary);	

	sf_floatwrite(output,n1,out);
    }

    elpf_close();

    exit (0);
}

/* 	$Id$	 */
