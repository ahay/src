/* Non-local smoothing. */
#include <rsf.h>

int main (int argc, char *argv[])
{
    int n1,n2, i1,i2, is, ns;
    float *trace, *trace2, ax, ay, t;
    sf_file inp, out;

    /* initialize */
    sf_init(argc,argv);

    /* set input and output files */
    inp = sf_input("in");
    out = sf_output("out");

    /* get input dimensions */
    if (!sf_histint(inp,"n1",&n1)) 
	sf_error("No n1= in input");
    n2 = sf_leftsize(inp,1);
    
    /* get command-line parameters */
    if (!sf_getint("ns",&ns)) sf_error("Need ns=");
    /* spray radius */

    if (!sf_getfloat("ax",&ax)) sf_error("Need ax=");
    /* exponential weight for the coordinate distance */

    trace = sf_floatalloc(n1);
    trace2 = sf_floatalloc(n1);

    /* loop over traces */
    for (i2=0; i2 < n2; i2++) {
	/* read input */
	sf_floatread(trace,n1,inp);

	/* loop over samples */
	for (i1=0; i1 < n1; i1++) {	    
	    t = 0.;

	    /* accumulate shifts */
	    for (is=-ns; is <= ns; is++) {
		if (i1+is >= 0 && i1+is < n1) {

		    /* !!!MODIFY THE NEXT LINE!!! */
		    t += trace[i1+is]*expf(-ax*is*is);
		} 
	    }

	    trace2[i1] = t;
	}

	/* write output */
	sf_floatwrite(trace2,n1,out);
    }

    /* clean up */
    sf_fileclose(inp);
    exit (0);
}


