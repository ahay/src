/* Clip the data. */

#include <rsf.h>

int main(int argc, char* argv[])
{
    int n, i, nclip;
    float clip, pclip, *data, *adata;
    sf_file in, out; /* Input and output files */

    /* Initialize RSF */
    sf_init(argc,argv);
    /* standard input */
    in  = sf_input("in");
    /* standard output */
    out = sf_output("out");

    /* check that the input is float */
    if (SF_FLOAT != sf_gettype(in)) 
	sf_error("Need float input");

    /* get data size */
    n = sf_filesize(in);

    /* parameter from the command line (i.e. pclip=90 ) */
    if (!sf_getfloat("pclip",&pclip)) sf_error("Need pclip=");

    /* allocate floating point arrays */
    data = sf_floatalloc (n);
    adata = sf_floatalloc (n);

    /* read entire dataset */
    sf_floatread(data,n,in);

    /* take absolute value */
    for (i=0; i < n; i++) {
	adata[i] = fabsf(data[i]);
    }

    if (pclip <= 0.0f) {
	/*
	  Be creative and find an optimal clip automatically! 
	*/
    } else {
	/* convert pclip to clip */
	nclip = SF_MAX(SF_MIN(n*pclip/100.0 + 0.5,n-1),0);
	clip = sf_quantile( nclip, n, adata);
    }

    /* loop over data samples */
    for (i=0; i < n; i++) {
	if      (data[i] >  clip) data[i]= clip;
	else if (data[i] < -clip) data[i]=-clip;
    }

    /* write the clipped data */
    sf_floatwrite(data,n,out);

    exit(0);
}
