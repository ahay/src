/* Clip the data.

Takes: < input.rsf > output.rsf 
*/
#include <rsf.h>

int main(int argc, char* argv[])
{
    int i, n, nbuf;
    float clip, *trace;
    sf_file in, out; /* Input and output files */

    /* Initialize RSF */
    sf_init(argc,argv);
    /* standard input */
    in = sf_input("in");
    /* standard output */
    out = sf_output("out");

    /* check that the input is float */
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    n = sf_filesize(in);

    /* parameter from the command line (i.e. clip=1.5 ) */
    if (!sf_getfloat("clip",&clip)) sf_error("Need clip=");

    /* allocate floating point array */
    nbuf = BUFSIZ/sizeof(float);
    trace = sf_floatalloc (nbuf);

    /* loop over traces */
    for (n = sf_filesize(in); n > 0; n -= nbuf) {
	if (nbuf > n) nbuf=n;

	/*read a trace */
	sf_floatread(trace,nbuf,in);

	/* loop over samples */
	for (i=0; i < n; i++) {
	    if      (trace[i] >  clip) trace[i]= clip;
	    else if (trace[i] < -clip) trace[i]=-clip;
	}
    
	/* write a trace */
	sf_floatwrite(trace,nbuf,out);
    }
    
    sf_close();
    exit(0);
}
