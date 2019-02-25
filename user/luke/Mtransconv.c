/* 
   program translates a 2D image then convolves it with arbitrary kernel
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
		/* variable translations file with same sampling as input, added ndim dimension */
		fixed = false;
		/* set translation file */
		_trans = sf_input("trans");
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
	
	/* convolution kernel boolean */
	bool ker;
	/* convolution kernel file */
	sf_file _ker;
	/* kernel dimension integers */
	int nk1, nk2;	
	/* kernel dimension array */
	int* NK ;
	/* kernel array */
	float* kernel;
	
	if ( NULL != sf_getstring("ker") ) {
	    /* convolution kernel file */
		/* set boolean */
		ker = true ;
	    /* set kernel */
	    _ker = sf_input("ker");
	    /* get kernel dimensions */
        if (!sf_histint  (_ker,"n1",&nk1))   sf_error("No n1= in kernel ");
        if (!sf_histint  (_ker,"n2",&nk2))   sf_error("No n2= in kernel");
	    /* allocate kernel dimension array */
	    NK = sf_intalloc(ndim);
   	    /* copy kernel size to array */
   	    NK[0] = nk1;
   	    NK[1] = nk2;
		/* allocate convolution kernel array */
		kernel = sf_floatalloc(conv_arraysize(NK,ndim));
		/* read convoluton kernel */
		sf_floatread(kernel,conv_arraysize(NK,ndim),_ker);
    } else {
		/* no convolution kernel, not convolving only translating */
    	ker = false ;
    }
	
	/* declare input array */
	float* arrayin  = sf_floatalloc(conv_arraysize(N,ndim));
	/* and output array */
	float* arrayout = sf_floatalloc(conv_arraysize(N,ndim));
	/* read from file */
	sf_floatread(arrayin,conv_arraysize(N,ndim),_in);
	
	if ( ker ){
		/* perform convolution */
		if (fixed){
			/* fixed translation convolution */
			arrayout =  conv_convolve_ker_translate( arrayin, X, N, D, O, kernel, NK, ndim, adj);
			/* free X */
			free (X);
		} else {
			/* variable translation convolution */
			arrayout = conv_convolve_ker_var_translate_omp( arrayin, trans, N, D, O, kernel, NK, ndim, adj);
			/* free translation array */
			free ( trans);
			/* close input file for translation */
			sf_fileclose(_trans);
		}
		/* free kernel array */
		free (kernel);
		/* and its dimensions */
		free (NK);
		/* and close the file */
		sf_fileclose(_ker);
	}else{
	    /* translate input, no convolution */
	    if (fixed){
		    /* translate by fixed coordinates */
		    arrayout = conv_translate_wrap(arrayin, X, N, D, O, ndim, adj);
			/* free X array */
			free ( X );
	    } else {
		    /* translate by variable coordinates */
		    arrayout = conv_var_translate_wrap(arrayin, trans, N, D, O, ndim, adj);
		    /* free translation array */
		    free(trans);
			/* close input file for translation */
			sf_fileclose(_trans);
	    }
    }
	/* write translated array to file */
	sf_floatwrite(arrayout,conv_arraysize(N,ndim),_out);
	/* close output file */
	sf_fileclose(_out);
	/* close input file */
	sf_fileclose(_in);
	/* free the substantial arrays*/
	free ( arrayin );
	free ( arrayout);
	/* and the smaller ones */
	free (N);
	free (D);
	free (O);

    /* exit program */
    exit (0);
}
