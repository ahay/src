/* Local signal and noise separation (N-dimensional).

Takes: < data.rsf sfilt=sigpef.rsf nfilt=noisepef.rsf > signal.rsf

Signal and noise separation by inversion (super-deconvolution).
Uses the helix and patching technologies.
*/

#include <rsf.h>

#include "helix.h"
#include "tent.h"
#include "patching.h"
#include "signoi.h"

int main(int argc, char* argv[])
{
    int n123, n1, i, ik, dim, nk, nf, sf, niter, nw;
    int n[SF_MAX_DIM], w[SF_MAX_DIM], k[SF_MAX_DIM];
    int sa[SF_MAX_DIM], na[SF_MAX_DIM], sc[SF_MAX_DIM], nc[SF_MAX_DIM];
    float *data, *wind, *sign, eps, di, dabs;
    char varname[6], *lagfile;
    filter saa, naa, sbb, nbb, scc, ncc;
    sf_file dat, signal, spef, npef, slag, nlag;

    sf_init (argc,argv);
    dat = sf_input("in");
    signal = sf_output("out");

    spef = sf_input("sfilt");
    npef = sf_input("nfilt");

    n123 = sf_filesize(dat);
    if (!sf_histint(spef,"dim",&dim)) sf_error("No dim= in sfilt");

    n1 = 1;
    for (i=0; i < dim; i++) {
	sprintf(varname,"n%d",i+1);
	if (!sf_histint(dat,varname,n+i)) 
	    sf_error("No %s= in input",varname);
	n1 *= n[i];
    }

    if (!sf_histints(spef,"w",w,dim)) sf_error("No w= in sfilt");
    if (!sf_histints(spef,"k",k,dim)) sf_error("No k= in sfilt");

    if (!sf_histints(spef,"a",sa,dim)) sf_error("No a= in sfilt");
    if (!sf_histints(npef,"a",na,dim)) sf_error("No a= in nfilt");

    if (!sf_histints(spef,"center",sc,dim)) sf_error("No center= in sfilt");
    if (!sf_histints(npef,"center",nc,dim)) sf_error("No center= in nfilt");

    nk=nw=1;
    for (i=0; i < dim; i++) {
	nw *= w[i];
	nk *= k[i];
    }

    if (!sf_histint(spef,"n1",&sf)) sf_error("No n1= in sfilt");
    if (!sf_histint(npef,"n1",&nf)) sf_error("No n1= in nfilt");

    sbb = allocatehelix(sf);
    nbb = allocatehelix(nf);

    if (NULL == (lagfile = sf_histstring(spef,"lag")) &&
	NULL == (lagfile = sf_getstring("slag"))) 
	sf_error("Need slag=");
    slag = sf_input(lagfile);
    if (NULL == (lagfile = sf_histstring(npef,"lag")) &&
	NULL == (lagfile = sf_getstring("nlag"))) 
	sf_error("Need nlag=");
    nlag = sf_input(lagfile);

    sf_intread(sbb->lag,sf,slag);
    sf_intread(nbb->lag,nf,nlag);

    if (!sf_getfloat("eps",&eps)) sf_error("Need eps=");
    /* regularization parameter */
    if (!sf_getint("niter",&niter)) niter=20;
    /* number of iterations */

    data = sf_floatalloc(n123);
    sign = sf_floatalloc(n123);

    sf_floatread(data,n123,dat);

    dabs = fabsf(data[0]);
    for (i=1; i < n123; i++) {
	di = fabsf(data[i]);
	if (di > dabs) dabs=di;
    }

    for (i=0; i < n123; i++) {
	data[i] /= dabs;
    }


    saa = (filter) sf_alloc(nk,sizeof(*saa));
    naa = (filter) sf_alloc(nk,sizeof(*naa));

    for (ik=0; ik < nk; ik++) {
	scc = saa+ik;
	ncc = naa+ik;
	scc->nh = sf;
	ncc->nh = nf;
	scc->flt = sf_floatalloc(sf);
	ncc->flt = sf_floatalloc(nf);
	scc->lag = sbb->lag;
	ncc->lag = nbb->lag;
	scc->mis = NULL;
	ncc->mis = NULL;
    }

    wind = sf_floatalloc(nw);
    tent (dim, w, (sc < nc)? sc: nc, (sa < na)? sa: na, wind);
 
    for (i=0; i < n123-n1+1; i += n1) {
	signoi_init (naa, saa, niter, nw, eps);
	for (ik=0; ik < nk; nk++) {
	    sf_floatread((naa+ik)->flt,nf,npef);
	    sf_floatread((saa+ik)->flt,sf,spef);
	}
	patching (signoi_lop, data+i, sign+i, dim, k, n, w, wind);
    }

    sf_floatwrite (sign,n123,signal);

    sf_close();
    exit(0);
}
