/* Time power gain 

Takes: < in.rsf > out.rsf

Does not estimate tpow automatically.
*/
#include <math.h>

#include <rsf.h>

int main(int argc, char *argv[])
{
    int nt, n2, i2, it;
    float tpow, *trace, *tgain=NULL, dt, t0;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&dt)) dt=1.;
    if (!sf_histfloat(in,"o1",&t0)) t0=0.;
    n2 = sf_leftsize(in,1);

    if (!sf_getfloat("tpow",&tpow)) tpow=2.;

    trace = sf_floatalloc(nt);

    if (0. != tpow) {
	tgain = sf_floatalloc(nt);
	for (it=0; it < nt; it++) {
	    tgain[it] = powf(t0+(it+1)*dt,tpow);
	}
    }

    for (i2=0; i2 < n2; i2++) {
	sf_read(trace,sizeof(float),nt,in);

	if (0. != tpow) {
	    for (it=0; it < nt; it++) {
		trace[it] *= tgain[it];
	    }
	}

	sf_write(trace,sizeof(float),nt,out);
    }

    sf_close();
    exit(0);
}
