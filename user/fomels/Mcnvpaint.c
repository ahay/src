/* Conection painting. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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
  int order, n1, n2, n12, i1, i2, i0, nw;
  char *label, *unit;
    float **u, ***conv, *trace;
    float o1, d1, o2, d2, eps, *time;
    sf_file cnv, out, seed;

    sf_init(argc, argv);
    cnv = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(cnv)) sf_error("Need float input");
    if (!sf_histint(cnv,"n1",&n1)) sf_error("Need n1= in input");
    
    if (!sf_histint(cnv,"n3",&n2)) sf_error("Need n3= in input");
    if (!sf_histfloat(cnv,"o3",&o2)) o2=0.; 
    if (!sf_histfloat(cnv,"d3",&d2)) d2=1.; 
    
    n12 = n1*n2;
    
    if (!sf_histint(cnv,"n2",&nw)) sf_error("Need n2= in input");

    sf_putint(out,"n2",n2);
    sf_putfloat(out,"d2",d2);
    sf_putfloat(out,"o2",o2);
    if (NULL != (label= sf_histstring(cnv,"label3"))) sf_putstring(out,"label2",label);
    if (NULL != (unit= sf_histstring(cnv,"unit3"))) sf_putstring(out,"unit2",unit);
    
    sf_putint(out,"n3",1);
    
    order = nw/2;

    if (NULL != sf_getstring("seed")) {
	seed = sf_input("seed");
	time = NULL;
    } else {
	seed = NULL;
	time = sf_floatalloc(n1);
	if (!sf_histfloat(cnv,"o1",&o1)) o1=0.; 
	if (!sf_histfloat(cnv,"d1",&d1)) d1=1.; 
	for (i1=0; i1 < n1; i1++) {
	    time[i1] = o1+i1*d1;
	}
    }

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */

    if (!sf_getint("i0",&i0)) i0=0;
    /* reference trace */
    
    predict_init (n1, n2, eps*eps, order, 1, false);

    u = sf_floatalloc2(n1,n2);
    conv = sf_floatalloc3(n1,nw,n2);
    trace = sf_floatalloc(n1);

    sf_floatread(conv[0][0],n12*nw,cnv);

    if (NULL != seed) {
      sf_floatread(trace,n1,seed);
    } else {
      for (i1=0; i1 < n1; i1++) {
	trace[i1] = time[i1];
      }
    }

    for (i1=0; i1 < n1; i1++) {
      u[i0][i1] = trace[i1];
    }
    for (i2=i0-1; i2 >= 0; i2--) {
      predict_step(false,false,trace,conv[i2]);
      for (i1=0; i1 < n1; i1++) {
	u[i2][i1] = trace[i1];
      }
    }
    for (i1=0; i1 < n1; i1++) {
      trace[i1] = u[i0][i1];
    }
    for (i2=i0+1; i2 < n2; i2++) {
      predict_step(false,true,trace,conv[i2-1]);
      for (i1=0; i1 < n1; i1++) {
	u[i2][i1] = trace[i1];
      }
    }
    sf_floatwrite(u[0],n1*n2,out);
  
    exit(0);
}
  

  
