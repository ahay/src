/* Signal and noise separation (N-dimensional).

Takes: < data.rsf sfilt=sigpef.rsf nfilt=noipef.rsf > signoi.rsf
*/

#include <rsf.h>

#include "helix.h"
#include "signoi.h"
#include "regrid.h"

int main(int argc, char* argv[])
{
    int niter, sa, na, j, dim, nx, n[SF_MAX_DIM], m[SF_MAX_DIM];
    float *dd, *ss, eps, na0, sa0;
    char varname[6], *lagfile;
    filter naa, saa;
    sf_file spef, npef, dat, signoi, slag, nlag;

    sf_init(argc,argv);
    dat = sf_input("in");
    signoi = sf_output("out");

    spef = sf_input("sfilt");
    npef = sf_input("nfilt");

    dim = sf_filedims(dat,n);
    
    if (!sf_histint(spef,"n1",&sa)) sf_error("No n1= in sfilt");
    if (!sf_histint(npef,"n1",&na)) sf_error("No n1= in nfilt");

    naa = allocatehelix(na);
    saa = allocatehelix(sa);

    if (NULL == (lagfile = sf_histstring(spef,"lag")) &&
	NULL == (lagfile = sf_getstring("slag"))) 
	sf_error("Need slag=");
    slag = sf_input(lagfile);
    if (!sf_histints(slag,"n",m,dim)) sf_error("No n= in %s",lagfile);
    sf_intread(saa->lag,sa,slag);
    regrid(dim,m,n,saa);
    sf_fileclose(slag);

    if (NULL == (lagfile = sf_histstring(npef,"lag")) &&
	NULL == (lagfile = sf_getstring("nlag"))) 
	sf_error("Need nlag=");
    nlag = sf_input(lagfile);
    if (!sf_histints(slag,"n",m,dim)) sf_error("No n= in %s",lagfile);
    sf_intread(naa->lag,na,nlag);
    regrid(dim,m,n,naa);
    sf_fileclose(nlag);

    if (!sf_histfloat(spef,"a0",&sa0)) sa0=1.;
    if (!sf_histfloat(npef,"a0",&na0)) na0=1.;

    if (!sf_getfloat("epsilon",&eps)) sf_error("Need eps=");
    /* regularization parameter */

    if (!sf_getint("niter",&niter)) niter=20;
    /* number of iterations */

    sprintf(varname,"n%d",dim+1);
    sf_putint(signoi,varname,2);

    nx=1;
    for(j=0; j < dim; j++) {
	nx *= n[j];
    }

    dd = sf_floatalloc(nx);
    ss = sf_floatalloc(nx);

    sf_floatread(dd,nx,dat);
    sf_floatread(saa->flt,sa,spef);
    for (j=0; j < sa; j++) {
	saa->flt[j] /= sa0;
    }

    sf_floatread(naa->flt,na,npef);
    for (j=0; j < na; j++) {
	naa->flt[j] /= na0;
    }

    signoi_init (naa, saa, niter, nx, eps);
    signoi_lop  (false,false,nx,nx,dd,ss);

    sf_floatwrite(ss,nx,signoi);

    for (j=0; j < nx; j++) {
	dd[j] -= ss[j];
    }

    sf_floatwrite(dd,nx,signoi);

    sf_close();
    exit (0);
}

