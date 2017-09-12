/* Streaming prediction filter in t-x domain. */
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

    int i1, i2, it, ix, n1, n2, n12, dim, na, i, nst;
    /* Define variables */
    
    int a[SF_MAX_DIM], n[SF_MAX_DIM];
    /* Define arrays for filter size a and data size n */
    
    float dd, da, dn, rn, lambda, lambda1, lambda2;
    /* Define variables */
    
    float *d, *r, *aa, *st;
    /* Define arrays for d, r, (bar)a_t, and (bar)a_x in equation 9 */

    sf_file in, out;
    /* Define file pointers of input and output  */
   
    sf_init(argc,argv);
    /* Initialize parameters for input and output */
    
    in = sf_input("in");
    out = sf_output("out");
    /* Initialize file pointers of input and output */
   
    dim = sf_filedims(in,n);
    /* Get data size and dimensions from input */
    if (2 < dim) dim = 2;
    
    if (!sf_getints("a",a,dim)) sf_error("Need a=");
    /* Get filter size from input, a0 is 2M+1, a1 is N in equation 3 */
    if (dim < 2) sf_error("Need at least two dimension");

    a[1]=a[1]*2;
    /* a1 is changed to 2N */
    
    n12 = 1;
    na = 1;
    
    for (i=0; i < dim; i++) {
	n12 *= n[i];
	na *= a[i];
    }
    
    n1=n[0];
    n2=n[1];
    nst=na*n1;
    
    if (!sf_getfloat("lambda1",&lambda1)) sf_error("Need lambda1=");
    /* Regularization in t direction, lambda_t in equations 1 and 5 */

    lambda1*=lambda1;
    /* Calculate square of lambda_t in equation 5 */

    if (!sf_getfloat("lambda2",&lambda2)) sf_error("Need lambda2=");
    /* Regularization in x direction, lambda_x in equations 1 and 5 */

    lambda2*=lambda2;
    /* Calculate square of lambda_x in equation 5 */

    lambda=lambda1+lambda2;
    /* Calculate square of lambda in equation 5 */
    
    d = sf_floatalloc(n12);
    /* Open space of array variable d(t,x) in equations 1, 8, 9 */
    r = sf_floatalloc(n12);
    /* Open space of array variable r(t,x) in equation 9 */
    aa = sf_floatalloc(na);
    /* Open space of array variable (bar)a_t in equation 5 */
    st = sf_floatalloc(nst);
    /* Open space of array variable (bar)a_x in equation 5 */
    
    for (i=0;i<na;i++){
	aa[i]=0.0f;
    }
    for (i=0;i<nst;i++){
	st[i]=0.0f;
    }
    
    sf_floatread(d,n12,in);
    /* Read data from input to array d(t,x) in equation 9 */
    
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    dd = 0.0f;
	    da = 0.0f;
	    i=0;
	    for (ix=-a[1]/2; ix < a[1]/2+1; ix++) {
		for (it=-a[0]/2; it < (a[0]+1)/2; it++) {
		    if(ix!=0){
			if(i2+ix<0 ||
			   i2+ix>=n2 ||
			   i1+it<0 ||
			   i1+it>=n1){
			    i++;
			    continue;
			} else{
			    dd += d[(i2+ix)*n1+i1+it]*
				d[(i2+ix)*n1+i1+it];
			    /* Variable dd is d^T d in equation 9 */
			    
			    da += d[(i2+ix)*n1+i1+it]*
				(lambda1*aa[i]+lambda2*st[i1*na+i])/lambda;
			    /* Variable da is d^T bar(a) in equation 9 */
			    i++;
			}
		    }
		}
	    }
	    /* Calculate d^T d and d^T bar(a) in equation 9 */
	    
	    dn=d[i2*n1+i1];
	    /* Variable dn is d(t,x) in equation 9 */

	    rn = (dn+da)/(lambda+dd);
	    r[i2*n1+i1] = lambda*rn;
	    /* Implement equation 9 */

	    i=0;
	    for (ix=-a[1]/2; ix < a[1]/2+1; ix++) {
		for (it=-a[0]/2; it < (a[0]+1)/2; it++) {
		    if(ix!=0){
			if(i2+ix<0 ||
			   i2+ix>=n2 ||
			   i1+it<0 ||
			   i1+it>=n1){
			    i++;
			    continue;
			} else{

			    aa[i] = (lambda1*aa[i]+lambda2*st[i1*na+i])/
				lambda-rn*d[(i2+ix)*n1+i1+it];
			    /* Implement equation 8 */
			    i++;
			}
		    }
		}
	    }
	    /* Calculate previous time-neighboring PF bar(a_t) in equation 5 */
	    for (i=0; i < na; i++) {
		st[i1*na+i]=aa[i];
	    }
	    /* Store previous space-neighboring PF column bar(a_x) in equation 5 */
	}
    }

    sf_floatwrite(r,n12,out);
    /* Output r(t,x) in equation 9 */

    exit(0);
}

