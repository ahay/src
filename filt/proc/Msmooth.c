/* Multi-dimensional smoothing.

Takes: < input.rsf > smooth.rsf rect1=1 rect2=1 ... 

rectN defines the size of the smoothing stencil in N-th dimension.
*/

#include <rsf.h>

#include "triangle.h"

static int first_index (int i, int j, int dim, const int *n, const int *s);

int main (int argc, char* argv[]) 
{
    int dim, dim1, i, j, n[SF_MAX_DIM], rect[SF_MAX_DIM], s[SF_MAX_DIM];
    int nrep, irep, n1, n2, i2, i0;
    bool diff[SF_MAX_DIM];
    char key[6];
    float* data;
    triangle tr;
    sf_file in, out;

    sf_init (argc, argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    dim = sf_filedims (in,n);
    dim1 = -1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	if (rect[i] > 1) dim1 = i;
	snprintf(key,6,"diff%d",i+1);
	if (!sf_getbool(key,diff+i)) diff[i]=false;
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

    if (!sf_getint("repeat",&nrep)) nrep=1;
    /* repeat filtering several times */

    for (i2=0; i2 < n2; i2++) {
	sf_read(data,sizeof(float),n1,in);

	for (i=0; i <= dim1; i++) {
	    if (rect[i] <= 1) continue;
	    tr = triangle_init (rect[i],n[i]);
	    for (j=0; j < n1/n[i]; j++) {
		i0 = first_index (i,j,dim1+1,n,s);
		for (irep=0; irep < nrep; irep++) {
		    smooth (tr,i0,s[i],diff[i],data);
		}
	    }
	    triangle_close(tr);
	}
	
	sf_write(data,sizeof(float),n1,out);
    }    

    exit (0);
}

static int first_index (int i, int j, int dim, const int *n, const int *s)
{
    int i0, n123, k, ii;

    n123 = 1;
    i0 = 0;
    for (k=0; k < dim; k++) {
	if (k == i) continue;
	ii = (j/n123)%n[k]; /* to cartesian */
	n123 *= n[k];	
	i0 += ii*s[k];      /* back to line */
    }

    return i0;
}
