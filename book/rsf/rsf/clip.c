#include <rsf.h>

int main(int argc, char* argv[])
{
    int n1, n2, i1, i2;
    float clip, *trace;
    sf_file in, out; 

    sf_init(argc,argv);
    in  = sf_input("in");
    out = sf_output("out");

    sf_histint(in,"n1",&n1); /* trace length */
    n2 = sf_leftsize(in,1);  /* number of traces */
    if (!sf_getfloat("clip",&clip)) sf_error("Need clip=");

    trace = sf_floatalloc (n1);
    for (i2=0; i2 < n2; i2++) {
	sf_floatread(trace,n1,in);
	for (i1=0; i1 < n1; i1++) {
	    if      (trace[i1] >  clip) trace[i1]= clip;
	    else if (trace[i1] < -clip) trace[i1]=-clip;
	}
	sf_floatwrite(trace,n1,out);
    }

    exit(0);
}
