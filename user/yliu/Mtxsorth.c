/* Streaming orthogonalize signal and noise in t-x domain. */
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

int main(int argc, char* argv[])
{
    int m[SF_MAX_DIM], i, it, ix, n1, n2, mdim, nd;
    /* Define variables */

    float gamma1, gamma2, remove;
    /* Define variables */

    float *s, *n, *w;
    /* Define arrays for s/s_0, n/n_0, w in equations 13 and 14 */

    sf_file in, noise, fnoi2, fsig2;
    /* Define file pointers of input and output  */

    sf_init(argc,argv);
    /* Initialize parameters for input and output */

    in = sf_input("in");
    noise = sf_input("noise");
    fnoi2 = sf_output("out");
    fsig2 = sf_output("sig2"); 
    /* Initialize file pointers of input and output */

    if (!sf_getfloat("gamma1",&gamma1)) sf_error("Need gamma1=");
    /* Regularization in t direction, gamma_t in equations 18 and 20 */

    gamma1 *= gamma1;
    /* Calculate square of gamma_t in equation 20 */

    if (!sf_getfloat("gamma2",&gamma2)) sf_error("Need gamma2=");
    /* Regularization in x direction, gamma_x in equations 18 and 20 */

    gamma2 *= gamma2;
    /* Calculate square of gamma_x in equation 20 */

    mdim = sf_filedims(in,m);
    /* Get data size and dimensions from input */

    nd = 1;
    for (i=0; i < mdim; i++) {
	nd *= m[i];
    }
    n1=m[0];
    n2=m[1];

    s = sf_floatalloc(nd);
    n = sf_floatalloc(nd);
    w = sf_floatalloc(nd);
    /* Open space of array variables for s, n, and w */

    sf_floatread(s,nd,in);
    sf_floatread(n,nd,noise);
    /* Read data arrays s_0 and n_0 from input */

    w[0]=0.0f;
    for(ix=0; ix < n2; ix++) {
	for(it=0; it < n1; it++) {
	    if(ix-1<0 || it-1<0) {
		w[ix*n1+it]=0.;
	    } else {
		w[ix*n1+it]=(s[ix*n1+it]*n[ix*n1+it]+
			     gamma1*w[ix*n1+it-1]+gamma2*w[(ix-1)*n1+it])/
		    (s[ix*n1+it]*s[ix*n1+it]+gamma1+gamma2);
		/* Implement equations 19 and 20 */
	    }
	}
    }
    
    for (i=0; i < nd; i++) {
	remove = w[i]*s[i];
	/* Implement w*s_0(t,x) in equations 13 and 14 */

	s[i] += remove;
        /* Implement equation 13 or 21 to calculate s(t,x) */

	n[i] -= remove;
        /* Implement equation 14 to calculate n(t,x) */
    }
    
    sf_floatwrite(n,nd,fnoi2);
    sf_floatwrite(s,nd,fsig2);
    /* Output s(t,x) and n(t,x) in equations 13 and 14 */
    
    exit(0);
}
