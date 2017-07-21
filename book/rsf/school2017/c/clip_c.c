/* Clip the data. */

#include <rsf.h>
#include <float.h>

int main(int argc, char* argv[])
{
    int n1, n2, i1, i2;
    float upper, lower;
    float *trace;
    /* Input and output files */
    sf_file in, out;

    /* Initialize RSF */
    sf_init(argc,argv);
    /* standard input */
    in  = sf_input("in");
    /* standard output */
    out = sf_output("out");

    /* check that the input is float */
    if (SF_FLOAT != sf_gettype(in))
	sf_error("Need float input");

    /* n1 is the fastest dimension (trace length) */
    if (!sf_histint(in,"n1",&n1)) 
	sf_error("No n1= in input");
    /* leftsize gets n2*n3*n4*... (the number of traces) */
    n2 = sf_leftsize(in,1);

    /* parameter from the command line (i.e. upper=1.) */
    if (!sf_getfloat("upper", &upper)) upper=+FLT_MAX;
    if (!sf_getfloat("lower", &lower)) lower=-FLT_MAX;

    /* allocate floating point array */
    trace = sf_floatalloc (n1);

    /* loop over traces */
    for (i2=0; i2 < n2; i2++) {

	/* read a trace */
	sf_floatread(trace,n1,in);

	/* loop over samples */
	for (i1=0; i1 < n1; i1++) {
	    if      (trace[i1] > upper) trace[i1]= upper;
	    else if (trace[i1] < lower) trace[i1]= lower;
	}

	/* write a trace */
	sf_floatwrite(trace,n1,out);
    }

    exit(0);
}
