/* 
   program calculates the alignment errors 
   between a reference and matching trace 
   given a maximum shift (maxshift=) and an error exponent (exp=)
*/

#include <rsf.h>
#include "dtw.h"

int main (int argc, char* argv[])
{
	/* declare files */
    sf_file _in, _ref, _out;
	/* initialize rsf */
    sf_init (argc,argv);
	/* get files from inputs */
    _in = sf_input("in");
    if ( NULL == sf_getstring("ref") ) { sf_error("Need reference trace ref=") ;}

    _ref  = sf_input ("ref");
	/* reference trace */ 
    _out = sf_output("out");
	int n1, maxshift;
	float d1, o1, ex;
    /* Get sampling info */
    if (!sf_histint  (_in,"n1",&n1))   sf_error("No n1=");
    if (!sf_histfloat(_in,"d1",&d1))   sf_error("No d1=");
	if (!sf_histfloat(_in,"o1",&o1))   sf_error("No o1=");

    if (!sf_getint("maxshift",&maxshift))   sf_error("Need maxshift=");
    /* maximum shift */
    if (!sf_getfloat("exp",&ex))   ex = 2;
    /* error exponent (g-f)^exp */
	if (ex < 1){ 
		sf_warning("Exponent doesn't define a norm, changing to exp=2");
		ex = 2;
	}
    /* allocate input matching trace array */
    float* match  = sf_floatalloc(n1);
	/* and read from file */
	sf_floatread(match,n1,_in);
	/* allocate secondary input reference trace array*/
	float* ref    = sf_floatalloc(n1);
	/* and read from file */
    sf_floatread( ref,n1,_ref);
	/* decleare mismatch array */
	float* mismatch = sf_floatalloc(n1*(2*maxshift+1));
    /* determine the best shifts using dynamic warping */
	dtw_alignment_errors( match, ref, mismatch, n1, maxshift, ex);
	/* spread values over nulls */
	dtw_spread_over_nulls( mismatch, n1, maxshift);
	/* set output sampling parameters */
	sf_putint   (_out,"n1",2*maxshift+1);
	sf_putfloat (_out,"d1",1);
	sf_putfloat (_out,"o1",-1.*maxshift);
	sf_putstring(_out,"label1", "lag");
	sf_putstring(_out, "unit1",  "l");
	sf_putint   (_out,"n2",n1);
	sf_putfloat (_out,"d2",d1);
	sf_putfloat (_out,"o2",o1);
	sf_putstring(_out,"label2", "sample index");
	sf_putstring(_out, "unit2",  "i");
	/* write the alignment errors array to output */
    sf_floatwrite( mismatch,n1*(2*maxshift+1),_out);
	/* deallocate arrays */
	free (mismatch);
	free (match   );
	free (ref     );
    /* exit program */
    exit (0);
}
