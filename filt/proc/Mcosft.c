/* Multi-dimensional cosine transform.

Takes: < input.rsf > transf.rsf sign1=0 sign2=0 ... 

signN defines transform along N-th dimension. 
signN=+1: forward transform.
signN=-1: inverse transform.

The input and output are real and have the same dimensions. 
Pad the data if you need to suppress wrap-around effects.
*/

#include <rsf.h>

#include "cosft.h"

int main (int argc, char* argv[]) 
{
    int dim, dim1, i, j, n[SF_MAX_DIM], sign[SF_MAX_DIM], s[SF_MAX_DIM];
    int n1, n2, i2, i0;
    char key[6], key2[6];
    float *data, o[SF_MAX_DIM], d[SF_MAX_DIM];
    sf_file in, out;

    sf_init (argc, argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    dim = sf_filedims (in,n);
    dim1 = -1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"sign%d",i+1);
	if (!sf_getint(key,sign+i)) sign[i]=0;
	if (sign[i]) {
	    dim1 = i;
	    snprintf(key,3,"o%d",i+1);
	    snprintf(key2,3,"u%d",i+1);
	    if (sign[i] > 0) {
		if (!sf_histfloat(in,key,o+i)) o[i]=0.;
		sf_putfloat(out,key,0.);
		sf_putfloat(out,key2,o[i]);
	    } else {		
		if (!sf_histfloat(in,key2,o+i)) o[i]=0.;
		sf_putfloat(out,key,o[i]);
	    }
	    snprintf(key,3,"d%d",i+1);
	    if (!sf_histfloat(in,key,d+i)) d[i]=1.;
	    sf_putfloat(out,key,1./(sf_npfar(2*(n[i]-1))*d[i]));
	}
    }

    n1 = n2 = 1;
    for (i=0; i < dim; i++) {
	if (i <= dim1) {
	    s[i] = n1;
	    n1 *= n[i];
	} else {
	    n2 *= n[i];
	}
    }

    data = sf_floatalloc (n1);

    for (i2=0; i2 < n2; i2++) {
	sf_read(data,sizeof(float),n1,in);
	
	for (i=0; i <= dim1; i++) {
	    if (!sign[i]) continue;
	    cosft_init(n[i] /* ,o[i],d[i] */);

	    for (j=0; j < n1/n[i]; j++) {
		i0 = sf_first_index (i,j,dim1+1,n,s);
		if (sign[i] > 0) {
		    cosft_frw (data, i0, s[i]);
		} else {
		    cosft_inv (data, i0, s[i]);
		}
	    }
	    cosft_close();
	}
	
	sf_write(data,sizeof(float),n1,out);
    }    
    
    exit (0);
}

/* 	$Id: Mcosft.c,v 1.4 2003/12/04 05:13:21 fomels Exp $	 */

