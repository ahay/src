/* Multidimensional convolution and deconvolution by helix transform.

Takes: < input.rsf filt=filter.rsf > output.rsf
*/

#include <rsf.h>

#include "helix.h"
#include "helicon.h"
#include "polydiv.h"
#include "regrid.h"

int main(int argc, char* argv[])
{
    int i, ia, na, nx, dim, n[SF_MAX_DIM], m[SF_MAX_DIM];
    float a0, *pp, *qq;
    bool adj, inv;
    filter aa;
    char* lagfile;
    sf_file in, out, filt, lag;

    sf_init (argc,argv);
    in = sf_input("in");
    filt = sf_input("filt");
    out = sf_output("out");

    dim = sf_filedims (in,n);

    if (!sf_histint(filt,"n1",&na)) sf_error("No n1= in filt");
    aa = allocatehelix (na);

    if (!sf_histfloat(filt,"a0",&a0)) a0=1.;
    sf_floatread (aa->flt,na,filt);
    for( ia=0; ia < na; ia++) {
	aa->flt[ia] /= a0;
    }

    if (NULL != (lagfile = sf_getstring("lag")) /* file with filter lags */
	|| 
	NULL != (lagfile = sf_histstring(filt,"lag"))) {
	lag = sf_input(lagfile);

	sf_intread(aa->lag,na,lag);
    } else {
	lag = NULL;
	for( ia=0; ia < na; ia++) {
	    aa->lag[ia] = ia+1;
	}
    }

    sf_fileclose(filt);
    
    if (!sf_getints ("n",m,dim) && (NULL == lag ||
				    !sf_histints (lag,"n",m,dim))) {
	for (i=0; i < dim; i++) {
	    m[i] = n[i];
	}
    }
 
    if (NULL != lag) sf_fileclose(lag);

    regrid (dim, m, n, aa);

    if (!sf_getbool ("adj",&adj)) adj=false;
    /* if y, do adjoint operation */
    if (!sf_getbool ("div",&inv)) inv=false;
    /* if y, do inverse operation (deconvolution) */

    nx = 1;
    for( i=0; i < dim; i++) {
	nx *= n[i];
    }
  
    pp = sf_floatalloc (nx);
    qq = sf_floatalloc (nx);

    if (adj) {
	sf_floatread (qq,nx,in);
    } else {
	sf_floatread (pp,nx,in);
    }

    if (inv) {
	polydiv_init (nx, aa);
	polydiv_lop (adj,false,nx,nx,pp,qq);
	polydiv_close();
    } else {
	helicon_init (aa);
	helicon_lop (adj,false,nx,nx,pp,qq);
    }

    if (adj) {
	sf_floatwrite (pp,nx,out);
    } else {
	sf_floatwrite (qq,nx,out);
    }

    sf_close();
    exit (0);
}

/* 	$Id: Mhelicon.c,v 1.5 2004/04/19 21:51:46 fomels Exp $	 */

