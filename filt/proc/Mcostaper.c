/* Cosine taper around the borders (2-D).

Takes: < in.rsf > tapered.rsf
*/
#include <math.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int n1, n2, nw, i1, i2, iw;
    float *trace, *w, wi;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getint("nw",&nw)) nw=10;
    /* tapering length */
    if (nw > n1) nw=n1;
    if (nw > n2) nw=n2;

    trace = sf_floatalloc(n1);
    w = sf_floatalloc(nw);

    for (iw=0; iw < nw; iw++) {
	wi = sinf(0.5*SF_PI*(iw+1.)/(nw+1.));
	w[iw] = wi*wi;
    }

    for (i2=0; i2 < n2; i2++) {
	sf_read(trace,sizeof(float),n1,in);
	for (iw=0; iw < nw; iw++) {
	    trace[iw] *= w[iw];
	    trace[n1-1-iw] *= w[iw];
	}
	if (i2 < nw) {
	    for (i1=0; i1 < n1; i1++) {
		trace[i1] *= w[i2];
	    }
	}
	if (n2-1-i2 < nw) {
	    for (i1=0; i1 < n1; i1++) {
		trace[i1] *= w[n2-1-i2];
	    }
	}
	sf_write(trace,sizeof(float),n1,out);
    }

    sf_close();
    exit(0);
}

