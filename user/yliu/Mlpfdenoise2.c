/* 2D denoising using edge-preserving local polynomial fitting (ELPF). */
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
    int n1, n2, n3, n12, nfw, i1, i2, i3, j, rect, niter, m, *coor, xc, yc, nw;
    bool boundary, verb;
    
    float *input, *output, *temp, *dd, *data;
    sf_file in, out, dip;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");
    dip = sf_input("dip");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    n12 = n1*n2;
    /* get the trace length (n1) and the number of traces (n2) and n3*/
    
    if (!sf_getint("nfw",&nfw)) sf_error("nfw needs integer input");
    /* filter-window length (positive and odd integer) */

    if (!sf_getint("nw",&nw)) sf_error("nw needs integer input");
    /* data-window length (positive and odd integer) */

    if (nw%2 ==0) nw +=1;
    m = (nw-1)/2;

    if (!sf_getint("rect",&rect)) sf_error("rect needs integer input");
    /* local smoothing radius */

    if (!sf_getbool("boundary",&boundary)) boundary = true;
    /* if y, boundary is data, whereas zero*/

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    
    input = sf_floatalloc(n12);
    output = sf_floatalloc(n12);
    dd = sf_floatalloc(n12);
    coor = sf_intalloc(nw*2);
    data = sf_floatalloc(nw);
    temp = sf_floatalloc(nw);

    elpf_init(nw,nfw,rect,false);
   
    for(i3=0; i3 < n3; i3++) {
	if (verb) sf_warning("3rd axis: %d of %d",i3+1,n3);
	sf_floatread(input,n12,in);
	sf_floatread(dd,n12,dip);

	for(i2=0; i2 < n2; i2++) {
	    if (verb && (10*i2)%n2==0) sf_warning("2nd axis: %d of %d",i2+1,n2);
	    for(i1=0; i1 < n1; i1++) {
		coor[0*nw+m] = i1;
		coor[1*nw+m] = i2;
		data[m] = input[coor[1*nw+m]*n1+coor[0*nw+m]];
		for(j=1; j < m+1; j++) {
		    xc = (int) coor[0*nw+j+m-1]+dd[coor[1*nw+j+m-1]*n1+coor[0*nw+j+m-1]]+0.5;
		    yc = coor[1*nw+j+m-1]+1;
		    if (xc < 0 || xc >= n1 || yc < 0 || yc >= n2) {
			coor[1*nw+j+m] = coor[1*nw+j+m-1];
			coor[0*nw+j+m] = coor[0*nw+j+m-1];
			if (boundary) {
			    data[j+m] = input[coor[1*nw+j+m]*n1+coor[0*nw+j+m]]; 
			} else {
			    data[j+m] = 0.;
			}
		    } else {
			coor[1*nw+j+m] = yc;
			coor[0*nw+j+m] = xc;
			data[j+m] = input[coor[1*nw+j+m]*n1+coor[0*nw+j+m]];
		    }
		}
		for(j=-1; j > (-1*(m+1)); j--) {
		    xc = (int) coor[0*nw+j+m+1]-dd[coor[1*nw+j+m+1]*n1+coor[0*nw+j+m+1]]+0.5;
		    yc = coor[1*nw+j+m+1]-1;
		    if (xc < 0 || xc >= n1 || yc < 0 || yc >= n2) {
			coor[1*nw+j+m] = coor[1*nw+j+m+1];
			coor[0*nw+j+m] = coor[0*nw+j+m+1];
			if (boundary) {
			    data[j+m] = input[coor[1*nw+j+m]*n1+coor[0*nw+j+m]]; 
			} else {
			    data[j+m] = 0.;
			}
		    } else {
			coor[1*nw+j+m] = yc;
			coor[0*nw+j+m] = xc;
			data[j+m] = input[coor[1*nw+j+m]*n1+coor[0*nw+j+m]];
		    }
		}
		output[i2*n1+i1] = elpf (data,temp,nw,nfw,niter,false);
	    }
	}
	sf_floatwrite(output,n12,out);
    }

    elpf_close();

    exit (0);
}

/* 	$Id$	 */
