#include <rsf.h>

#include "mshelicon.h"
#include "mshelix.h"
#include "helicon.h"

static msfilter aa;

/*
  Helicon
  -------
  Helical convolution. 
  Initialized with the filter. */
void mshelicon_init( msfilter bb) {
    aa = bb;
}

void mshelicon_lop( bool adj, bool add, int nx, int ny, float* xx, float*yy) {
    int is;
    
    sf_adjnull( adj, add, nx, ny, xx, yy);
    
    for (is=0; is < aa->ns; is++) {
	onescale(is,aa);
	helicon_init(aa->one);
	helicon_lop(adj,true,nx,nx,xx,yy+is*nx);
    }
}

/* 	$Id: mshelicon.c,v 1.1 2004/06/11 10:51:33 fomels Exp $	 */
