/* Find PEF on aliased traces.

Takes: < data.rsf > pef.rsf
*/

#include <rsf.h>

#include "helix.h"
#include "lace.h"
#include "printfilter.h"
#include "compress.h"

int main(int argc, char* argv[])
{
    int jump, dim, i, np;
    int n[SF_MAX_DIM], a[SF_MAX_DIM], center[SF_MAX_DIM], gap[SF_MAX_DIM];
    float *pp;
    char varname[6], *lagfile;
    filter aa;
    sf_file dat, pef, lag;

    sf_init(argc,argv);
    dat = sf_input("in");
    pef = sf_output("out");

    if (NULL == (lagfile = sf_getstring("lag"))) sf_error("Need lag=");
    /* output file for filter lags */

    lag = sf_output(lagfile);
    sf_settype(lag,SF_INT);

    sf_putstring(pef,"lag",lagfile);

    dim = sf_filedims(dat,n);

    if (!sf_getint("jump",&jump)) jump=2;

    if (!sf_getints("a",a,dim)) sf_error("Need a=");

    if (!sf_getints("center",center,dim)) {
	for (i=0; i < dim; i++) {
	    center[i] = (i+1 < dim && a[i+1] > 1)? a[i]/2: 0;
	}
    }

    if (!sf_getints("gap",gap,dim)) {
	for (i=0; i < dim; i++) {
	    gap[i] = 0;
	}
    }

    np = 1;
    for (i=0; i < dim; i++) {
	np *= n[i];
    }

    pp = sf_floatalloc(np);
    sf_floatread(pp,np,dat);

    aa = lace_pef (dim, pp, jump, np, n, center, gap, a);
    aa = compress(aa, 1.e-6);

    sf_putint(pef,"n1",aa->nh);
    sf_putint(lag,"n1",aa->nh);

    for (i=1; i < dim; i++) {
	sprintf(varname,"n%d",i+1);
	sf_putint(pef,varname,1);
	sf_putint(lag,varname,1);
	n[i] *= jump;
    }

    sf_putints (lag,"n",n,dim);
    sf_intwrite(aa->lag,aa->nh,lag);
    sf_floatwrite(aa->flt,aa->nh,pef);

    print (dim, n, center, a, aa);

    sf_close();
    exit(0);
}
