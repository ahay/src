/* Dip estimation from convection */
/*
  Copyright (C) 2025 University of Texas at Austin
  
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
#include <math.h>

#include <rsf.h>
#include "predict.h"

int main (int argc, char* argv[])
{
    int order, n1, n2, i1, i2, nw;
    float **conv, *trace, eps;
    sf_file cnv, dip;

    sf_init(argc, argv);
    cnv = sf_input("in");
    dip = sf_output("out");

    if (SF_FLOAT != sf_gettype(cnv)) sf_error("Need float input");
    if (!sf_histint(cnv,"n1",&n1)) sf_error("Need n1= in input");
    if (!sf_histint(cnv,"n2",&nw)) sf_error("Need n2= in input");    
    if (!sf_histint(cnv,"n3",&n2)) sf_error("Need n3= in input");

    order = nw/2;

    sf_unshiftdim(cnv,dip,2);

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */
    
    predict_init (n1, n2, eps*eps, order, 1, false);

    conv = sf_floatalloc2(n1,nw);
    trace = sf_floatalloc(n1);
 
    for (i2=0; i2 < n2; i2++) {	
	for (i1=0; i1 < n1; i1++) {
	    trace[i1] = i1;
	}
	sf_floatread(conv[0],n1*nw,cnv);
	predict_step(false,true,trace,conv);
	for (i1=0; i1 < n1; i1++) {
	    trace[i1] = i1 - trace[i1];
	}
	sf_floatwrite(trace,n1,dip);
    }
  
    exit(0);
}
  

  
