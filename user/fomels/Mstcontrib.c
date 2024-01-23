/* Contribution weighting using streaming attributes. */
/*
  Copyright (C) 2023 University of Texas at Austin
  
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

int main (int argc, char* argv[])
{
    int i1, n1, i2, n2;
    float w, eps, w0;
    float *stack, *gather, *weight;
    sf_file gath, stck, wght, outp;

    sf_init (argc,argv);
    gath = sf_input("in");
    stck = sf_input("stack");
    wght = sf_output("weight");
    outp = sf_output("out");

    if (!sf_histint(gath,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(gath,1);

    stack = sf_floatalloc(n1);
    gather = sf_floatalloc(n1);
    weight = sf_floatalloc(n1);

    if (!sf_getfloat("eps",&eps)) eps=1.0f;
    /* regularization parameter */
    eps *= eps;

    if (!sf_getfloat("w0",&w0)) w0=0.0f;
    /* initial weight */
    
    sf_floatread(stack,n1,stck);
    
    for (i2=0; i2 < n2; i2++) {
	w = w0;
	sf_floatread(gather,n1,gath);
	
	for (i1=0; i1 < n1; i1++) {
	    w = (gather[i1]*stack[i1]+eps*w)/(gather[i1]*gather[i1]+eps);
	    weight[i1]=w;
	    gather[i1]=w*gather[i1];
	}

	sf_floatwrite(weight,n1,wght);
	sf_floatwrite(gather,n1,outp);
    }

    exit(0);
}
