#include "triangle.h"
#include "trianglen.h"

#include <rsf.h>

static int *n, s[SF_MAX_DIM], nd, dim;
static triangle *tr;
static float *tmp;

void trianglen_init (int ndim, int *nbox, int *ndat)
{
    int i;

    n = ndat;
    dim = ndim;

    tr = (triangle*) sf_alloc(dim,sizeof(triangle));

    nd = 1;
    for (i=0; i < dim; i++) {
	tr[i] = (nbox[i] > 1)? triangle_init (nbox[i],ndat[i]): NULL;
	s[i] = nd;
	nd *= ndat[i];
    }
    tmp = sf_floatalloc (nd);
}

void trianglen_lop (bool adj, bool add, int nx, int ny, float* x, float* y)
{
    int i, j, i0;

    if (nx != ny || nx != nd) 
	sf_error("%s: Wrong data dimensions: nx=%d, ny=%d, nd=%d",
		 __FILE__,nx,ny,nd);

    sf_adjnull (adj,add,nx,ny,x,y);
  
    if (adj) {
	for (i=0; i < nd; i++) {
	    tmp[i] = y[i];
	}
    } else {
	for (i=0; i < nd; i++) {
	    tmp[i] = x[i];
	}
    }

  
    for (i=0; i < dim; i++) {
	if (NULL != tr[i]) {
	    for (j=0; j < nd/n[i]; j++) {
		i0 = sf_first_index (i,j,dim+1,n,s);
		smooth2 (tr[i], i0, s[i], false, tmp);
	    }
	}
    }
	
    if (adj) {
	for (i=0; i < nd; i++) {
	    x[i] += tmp[i];
	}
    } else {
	for (i=0; i < nd; i++) {
	    y[i] += tmp[i];
	}
    }    
}

void trianglen_close(void)
{
    int i;

    free (tmp);

    for (i=0; i < dim; i++) {
	if (NULL != tr[i]) triangle_close (tr[i]);
    }

    free(tr);
}

/* 	$Id: trianglen.c,v 1.1 2004/04/05 14:38:29 fomels Exp $	 */
