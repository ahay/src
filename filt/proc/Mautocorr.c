/* Autocorrelation for helix filters. */

#include <rsf.h>

#include "helix.h"
#include "autocorr.h"

int main(int argc, char* argv[])
{
    int i, na, ns;
    float s0, a0;
    char* lagfile;
    filter ss, aa;
    sf_file in, out, lag0, lag;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&na)) sf_error("No n1= in input");
    aa = allocatehelix (na);

    if (NULL == (lagfile = sf_histstring(in,"lag"))) {
	if (NULL == (lagfile = sf_getstring("lag"))) {
	    /* optional input file with filter lags */
	    for (i=0; i < na; i++) {
		aa->lag[i]=i+1;
	    }
	    lag0 = NULL;
	} else {
	    lag0 = sf_input("lag");
	}
    } else {
	lag0 = sf_input(lagfile);
    }

    if (NULL != lag0) {
	if (SF_INT != sf_gettype(lag0)) 
	    sf_error("Need int data in lag file '%s'",lagfile);

	sf_intread(aa->lag,na,lag0);
    }

    if (!sf_histfloat(in,"a0",&a0)) a0=1.;
    sf_floatread (aa->flt,na,in);

    ss = autocorr (aa, a0, &s0, 1.e-6);
    ns = ss->nh;

    sf_putfloat (out,"a0",s0*0.5);
    sf_putint (out,"n1",ns);

    lag = sf_output("lagout");
    sf_putint(lag,"n1",ns);
    sf_settype(lag,SF_INT);
    sf_fileflush(lag,lag0);

    sf_intwrite(ss->lag,ns,lag);
    sf_fileclose(lag);

    if (NULL != (lagfile = sf_getstring("lagout"))) 
	sf_putstring(out,"lag",lagfile);

    sf_floatwrite (ss->flt,ns,out);
  
    sf_close();
    exit (0);
}
