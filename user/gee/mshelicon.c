/* Multi-scale helical convolution */

#include <rsf.h>
/*^*/

#include "mshelicon.h"

#include "mshelix.h"
/*^*/

static msfilter aa;

void mshelicon_init( msfilter bb) 
/*< initialize with the filter >*/
{
    aa = bb;
}

void mshelicon_lop( bool adj, bool add, int nx, int ny, float* xx, float*yy) 
/*< linear operator >*/
{
    int is;
    
    sf_adjnull( adj, add, nx, ny, xx, yy);
    
    for (is=0; is < aa->ns; is++) {
	onescale(is,aa);
	sf_helicon_init(aa->one);
	sf_helicon_lop(adj,true,nx,nx,xx,yy+is*nx);
    }
}

/* 	$Id: mshelicon.c 2523 2007-02-02 16:45:29Z sfomel $	 */
