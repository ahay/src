/* Multi-dimensional smoothing with boxes.

Takes: < input.rsf > smooth.rsf rect1=1 rect2=1 ... 

rectN defines the size of the smoothing stencil in N-th dimension.
*/

#include <rsf.h>

#include "box.h"

int main (int argc, char* argv[]) 
{
    int dim, dim1, i, j, k, n[SF_MAX_DIM], rect[SF_MAX_DIM], s[SF_MAX_DIM];
    int nrep, irep, n1, n2, i2, i0, rmax;
    char key[6];
    float *data, *data2;
    sf_file in, out;

    sf_init (argc, argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    dim = sf_filedims (in,n);
    dim1 = -1;
    rmax = 1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	if (rect[i] > 1) {
	    dim1 = i;
	    if (rect[i]+n[i] > rmax) rmax = rect[i]+n[i];
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
    data2 = sf_floatalloc (rmax);

    if (!sf_getint("repeat",&nrep)) nrep=1;
    /* repeat filtering several times */

    for (i2=0; i2 < n2; i2++) {
	sf_read(data,sizeof(float),n1,in);

	for (i=0; i <= dim1; i++) {
	    if (rect[i] <= 1) continue;
	    box_init (rect[i],n[i],false);
	    for (j=0; j < n1/n[i]; j++) {
		i0 = sf_first_index (i,j,dim1+1,n,s);
		for (irep=0; irep < nrep; irep++) {
		    boxsmooth2 (i0,s[i],data,data2);
		    for (k=0; k < n[i]; k++) {
			data[i0+k*s[i]] = data2[k];
		    }
		}
	    }
	}
	
	sf_write(data,sizeof(float),n1,out);
    }    

    sf_close();
    exit (0);
}

/* 	$Id: Mbox.c,v 1.2 2004/03/22 05:43:24 fomels Exp $	 */
