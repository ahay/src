/* Cosine taper around the borders (2-D).

Takes: < in.rsf > tapered.rsf
*/
#include <math.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int n1, n2, nw1, nw2, i1, i2, iw;
    float *trace, *w1=NULL, *w2=NULL, wi;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getint("nw1",&nw1)) nw1=0;
    /* tapering length on first axis */
    if (nw1 > n1) nw1=n1;

    if (!sf_getint("nw2",&nw2)) nw2=0;
    /* tapering length on second axis */
    if (nw2 > n2) nw2=n2;

    trace = sf_floatalloc(n1);
    if (nw1 > 0) w1 = sf_floatalloc(nw1);
    if (nw2 > 0) w2 = sf_floatalloc(nw2);

    for (iw=0; iw < nw1; iw++) {
	wi = sinf(0.5*SF_PI*(iw+1.)/(nw1+1.));
	w1[iw] = wi*wi;
    }

    for (iw=0; iw < nw2; iw++) {
	wi = sinf(0.5*SF_PI*(iw+1.)/(nw2+1.));
	w2[iw] = wi*wi;
    }

    for (i2=0; i2 < n2; i2++) {
	sf_read(trace,sizeof(float),n1,in);
	for (iw=0; iw < nw1; iw++) {
	    wi = w1[iw];
	    trace[iw]      *= wi;
	    trace[n1-1-iw] *= wi;
	}
	if (i2 < nw2) {
	    wi = w2[i2];
	    for (i1=0; i1 < n1; i1++) {
		trace[i1] *= wi;
	    }
	}
	if (n2-1-i2 < nw2) {
	    wi = w2[n2-1-i2];
	    for (i1=0; i1 < n1; i1++) {
		trace[i1] *= wi;
	    }
	}
	sf_write(trace,sizeof(float),n1,out);
    }

    sf_close();
    exit(0);
}
