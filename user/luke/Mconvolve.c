/* 
   convolve input 2D image by kernel
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
	

	if ( NULL == sf_getstring("ker") ) sf_error( "Need ker=");
	/* convolution kernel file */
		
	/* convolution kernel */
	sf_file _ker = sf_input("ker");
	/* get kernel dimensions */
	int nk1, nk2;
    if (!sf_histint  (_ker,"n1",&nk1))   sf_error("No n1= in kernel ");
    if (!sf_histint  (_ker,"n2",&nk2))   sf_error("No n2= in kernel");

	int* NK = sf_intalloc(ndim);
	/* copy kernel size to array */
	NK[0] = nk1;
	NK[1] = nk2;
	
	/* declare input array */
	float* arrayin  = sf_floatalloc(conv_arraysize(N,ndim));
	/* and output array */
	float* arrayout = sf_floatalloc(conv_arraysize(N,ndim));
	/* read from file */
	sf_floatread(arrayin,conv_arraysize(N,ndim),_in);
	/* declare convolution kernel array */
	float* kernel = sf_floatalloc(conv_arraysize(NK,ndim));
	/* read convoluton kernel */
	sf_floatread(kernel,conv_arraysize(NK,ndim),_ker);
	/* convolve by kernel */
	arrayout = conv_convolve_ker( arrayin, N, kernel, NK, ndim, adj);
	/* write translated array to file */
	sf_floatwrite(arrayout,conv_arraysize(N,ndim),_out);
	/* close the output file */
	sf_fileclose(_out);
	/* close the input files */
	sf_fileclose(_in);
	sf_fileclose(_ker);
	/* free the substantial arrays*/
	free ( arrayin );
	free ( arrayout);
	free ( kernel );
	free (O);
	free (D);
	free (N);
	free (NK);
    /* exit program */
    exit (0);
}
