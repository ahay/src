/* 
   program translates a 2D image
*/

#include <rsf.h>
#include "convtools.h"

int main (int argc, char* argv[])
{
	/* how many dimensions? */
	int ndim = 2;
	/* declare files */
    sf_file _in, _out;
	/* initialize rsf */
    sf_init (argc,argv);
	/* get files from inputs */
    _in = sf_input("in");
	/* output */
    _out = sf_output("out");
	int n1, n2;
	float d1, d2, o1, o2;
	/* adjoint flag */
	bool adj;
    /* Get sampling info */
    if (!sf_histint  (_in,"n1",&n1))   sf_error("No n1=");
    if (!sf_histfloat(_in,"d1",&d1))   sf_error("No d1=");
	if (!sf_histfloat(_in,"o1",&o1))   sf_error("No o1=");
    if (!sf_histint  (_in,"n2",&n2))   sf_error("No n2=");
    if (!sf_histfloat(_in,"d2",&d2))   sf_error("No d2=");
	if (!sf_histfloat(_in,"o2",&o2))   sf_error("No o2=");
	
	if (!sf_getbool("adj",&adj)) adj=false;
	/*if y reverse translation, if n, translation */
	
	/* make N array */
	int* N = sf_intalloc(ndim);
	N[0] = n1;
	N[1] = n2;
	
	/* make D array */
	float* D = sf_floatalloc(ndim);
	D[0] = d1;
	D[1] = d2;
	
	/* make O array */
	float* O = sf_floatalloc(ndim);
	O[0] = o1;
	O[1] = o2;
	
	/* variable translation */
	sf_file _trans;
	/* translations array */
	float* trans;
	/* fixed translations */
	float x1, x2;
	/*  translation array */
	float* X;

	/* translations boolean */
	bool fixed;
	if ( NULL != sf_getstring("trans") ) { 
		fixed = false;
		_trans = sf_input("trans");
		/* variable translations file with same sampling as input, added ndim dimension */
		trans = sf_floatalloc( ((long)ndim) * conv_arraysize(N,ndim));
		/* read array */
		sf_floatread(trans, ((long)ndim) * conv_arraysize(N,ndim), _trans);
	} else {
		fixed = true;
		/* fixed translation */
	    if (!sf_getfloat("x1",&x1))   x1 = 0.;
		/* fixed translation in first dimension */
		if (!sf_getfloat("x2",&x2))   x2 = 0.;
		/* fixed translation in second dimension */
		
		/* allocate fixed translation array */
		X = sf_floatalloc(ndim);
		X[0] = x1;
		X[1] = x2;
	}
	
	/* declare input array */
	float* arrayin  = sf_floatalloc(conv_arraysize(N,ndim));
	/* and output array */
	float* arrayout = sf_floatalloc(conv_arraysize(N,ndim));
	/* read from file */
	sf_floatread(arrayin,conv_arraysize(N,ndim),_in);
	/* translate input */
	if (fixed){
		/* translate by fixed coordinates */
		arrayout = conv_translate_wrap(arrayin, X, N, D, O, ndim, adj);
	} else {
		/* translate by variable coordinates */
		arrayout = conv_var_translate_wrap(arrayin, trans, N, D, O, ndim, adj);
		/* free translation array */
		free(trans);
	}
	/* write translated array to file */
	sf_floatwrite(arrayout,conv_arraysize(N,ndim),_out);
	
	/* free the substantial arrays*/
	free ( arrayin );
	free ( arrayout);

    /* exit program */
    exit (0);
}
