/* 
   find kernel for convolution in ls sense, this is assuming 3 dimensional data and a 2d kernel
*/

#include <rsf.h>
#include "convtools.h"

int main (int argc, char* argv[])
{
	/* how many dimensions? */
	int ndim = 3;
	/* files */
	sf_file _in, _out, _match;
	/* initialize rsf */
    sf_init (argc,argv);
	/* adjoint flag */
	bool adj;
	if (!sf_getbool("adj",&adj)) adj=false;
	/*if n takes kernel and outputs function, if y takes function and outputs kernel */
	
	/* get files from inputs */
    _in = sf_input("in");
	/* output */
    _out = sf_output("out");
	
	if ( NULL == sf_getstring("match") ) sf_error( "Need match=");
	_match = sf_input("match");
	/* matching file we are mapping to, same dimension as function */
	
	/* sampling info */
	int n1, n2, n3;
	int nf1, nf2, nf3;
	int nk1, nk2;
	float d1, d2, d3;
	float o1, o2, o3;
	
    /* Get sampling info */
    if (!sf_histint  (_match,"n1",&n1))   sf_error("No n1= in match");
    if (!sf_histfloat(_match,"d1",&d1))   sf_error("No d1= in match");
	if (!sf_histfloat(_match,"o1",&o1))   sf_error("No o1= in match");
    if (!sf_histint  (_match,"n2",&n2))   sf_error("No n2= in match");
    if (!sf_histfloat(_match,"d2",&d2))   sf_error("No d2= in match");
	if (!sf_histfloat(_match,"o2",&o2))   sf_error("No o2= in match");
    if (!sf_histint  (_match,"n3",&n3))   sf_error("No n3= in match");
    if (!sf_histfloat(_match,"d3",&d3))   sf_error("No d3= in match");
	if (!sf_histfloat(_match,"o3",&o3))   sf_error("No o3= in match");
	
	/* make dimension arrays */
	/* make N array */
	int* N = sf_intalloc(ndim);
	N[0] = n1;
	N[1] = n2;
	N[2] = n3;
	
	/* make D array */
	float* D = sf_floatalloc(ndim);
	D[0] = d1;
	D[1] = d2;
	D[2] = d3;
	
	/* make O array */
	float* O = sf_floatalloc(ndim);
	O[0] = o1;
	O[1] = o2;
	O[2] = o3;
	
	/* allcoate match */
	float* match = sf_floatalloc(conv_arraysize(N,ndim-1));
	/* allocate function */
	float* function = sf_floatalloc(conv_arraysize(N,ndim-1));
	
	if (!adj){
		/* get kernel sampling from input file */
	    if (!sf_histint  (_in,"n1",&nk1))   sf_error("No n1= in input");
	    if (!sf_histint  (_in,"n2",&nk2))   sf_error("No n2= in input");

		
    	sf_putint   (_out,"n1",N[0]); 
	    sf_putfloat (_out,"d1",D[0]);
	    sf_putfloat (_out,"o1",O[0]);	
	    sf_putint   (_out,"n2",N[1]); 
	    sf_putfloat (_out,"d2",D[1]);
	    sf_putfloat (_out,"o2",O[1]);
	    sf_putint   (_out,"n3",N[2]);
	    sf_putfloat (_out,"d3",D[2]);
	    sf_putfloat (_out,"o3",O[2]);
		
	}else{
		/* user selected kernel info */
		
		if (!sf_getint("nk1",&nk1))   nk1 = 5;
		/* size of kernel in dimension 1*/
		
		if (!sf_getint("nk2",&nk2))   nk2 = 5;
		/* size of kernel in dimension 2 */
		
    	sf_putint   (_out,"n1",nk1); 
	    sf_putfloat (_out,"d1",1);
	    sf_putfloat (_out,"o1",-1*(nk1-1)/2);	
	    sf_putint   (_out,"n2",nk2); 
	    sf_putfloat (_out,"d2",1);
	    sf_putfloat (_out,"o2",-1*(nk2-1)/2);
	    sf_putint   (_out,"n3",1);
	    sf_putfloat (_out,"d3",1);
	    sf_putfloat (_out,"o3",0);
		
	    if (!sf_histint  (_in,"n1",&nf1))   sf_error("No n1= in input");
	    if (!sf_histint  (_in,"n2",&nf2))   sf_error("No n2= in input");
		if (!sf_histint  (_in,"n3",&nf3))   sf_error("No n3= in input");
		if (nf1 != n1 ) sf_error("Need n1 of input and match the same");
		if (nf2 != n2 ) sf_error("Need n2 of input and match the same");
		if (nf3 != n3 ) sf_error("Need n3 of input and match the same");		
	}
	
	int* Nk = sf_intalloc(ndim-1);
	Nk[0] = nk1;
	Nk[1] = nk2;
	
	/* allocate kernel */
	float* kernel = sf_floatalloc(conv_arraysize(Nk,ndim-1));
	
	/* read kernel if adjoint */
	if (!adj){
	    sf_floatread(kernel,conv_arraysize(Nk,ndim-1),_in);	
    }
	/* loop through 3rd dimension */
	int i3;
	for (i3 = 0 ; i3 < N[2] ; i3++){
		/* read function if not adjoint */
		if (adj) sf_floatread(function,conv_arraysize(N,ndim-1),_in);
		/* read match */
		sf_floatread(match,conv_arraysize(N,ndim-1),_match);
		conv_convolve_find_ker(match, N, kernel, function, Nk, ndim-1, adj)	;
		/* write function if appropriate */
		if (!adj) sf_floatwrite(function,conv_arraysize(N,ndim-1),_out);
		/* zero out array */
		conv_clear_array( function, conv_arraysize(N,ndim-1) );
	}
	/* write kernel */
	if (adj) sf_floatwrite(kernel,conv_arraysize(Nk,ndim-1),_out);
    /* exit program */
    exit (0);
}
