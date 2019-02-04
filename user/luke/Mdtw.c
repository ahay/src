/* 
   program calculates the shifts to warp a matching trace (input) 
   to a reference trace (ref=) of equal length and applies those shifts 
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
	int n1, maxshift, i;
	float d1, o1, ex, str;
    /* Get sampling info */
    if (!sf_histint  (_in,"n1",&n1))   sf_error("No n1=");
    if (!sf_histfloat(_in,"d1",&d1))   sf_error("No d1=");
	if (!sf_histfloat(_in,"o1",&o1))   sf_error("No o1=");
    /* maximum shift */
    /* Get smoothing radius */
    if (!sf_getint("maxshift",&maxshift))   sf_error("Need maxshift=");
    /* maximum shift */
    if (!sf_getfloat("exp",&ex))   ex = 2;
    /* error exponent (g-f)^exp */
	if (ex < 1){ 
		sf_warning("Exponent doesn't define a norm, changing to exp=2");
		ex = 2;
	}
    if (!sf_getfloat("strain",&str))   str = 1.0;
    /* maximum strain */
	if (str > 1){ 
		sf_warning("Program not set up for du/dt > 1, changing to str = 1");
		str = 1.0;
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
	/* declare array for forward error accumulation */
	float* accumulate_f = sf_floatalloc(n1*(2*maxshift+1));
	/* accumulate forward errors */
	dtw_accumulate_errors( mismatch, accumulate_f, n1, maxshift, str);
	/* declare array for backward error accumulation */
	float* accumulate_b = sf_floatalloc(n1*(2*maxshift+1));
	/* reverse mismatch */
	float* mismatch_r = sf_floatalloc(n1*(2*maxshift+1));
	dtw_reverse_array(mismatch, mismatch_r, n1, (2*maxshift+1));
	/* accumulate backward errors */
	dtw_accumulate_errors( mismatch_r, accumulate_b, n1, maxshift, str);
	free (   mismatch_r);
	/* flip them */
	float* accumulate_br = sf_floatalloc(n1*(2*maxshift+1));
	dtw_reverse_array(accumulate_b,accumulate_br, n1, (2*maxshift+1));
	free ( accumulate_b) ;
	/* sum the errors */
	float* accumulate = sf_floatalloc(n1*(2*maxshift+1));
	dtw_two_sided_smoothing_sum(accumulate_f, accumulate_br, mismatch, accumulate,n1*(2*maxshift+1));
	free ( accumulate_f) ;
	free ( accumulate_br);
	/* allocate shifts array */
	int* shifts = sf_intalloc(n1);
	/* backtrack to find integer shifts */
	dtw_backtrack( accumulate, mismatch, shifts, n1, maxshift, str);
	/* are we writing out the mismatch?*/
	sf_file _error;
    if ( NULL != sf_getstring("error") ) {
	/* misfit error */ 
	    _error  = sf_output ("error");
		/* set output parameters */
		dtw_set_file_params(_error, 2*maxshift+1, 1, -1*(float)maxshift,"lag","l", n1, d1, o1,"sample index","i");
	    /* write to file */
    	sf_floatwrite(mismatch,n1*(2*maxshift+1),_error);
    }
	/* free mismatch array */
	free (   mismatch  );
	
	/* are we writing out the accumulation errors? */
	sf_file _accum;
    if ( NULL != sf_getstring("accum") ) {
	/* accumulation errors from forward and backtracking */ 
	    _accum  = sf_output ("accum");
		/* set output parameters */
		dtw_set_file_params(_accum, 2*maxshift+1, 1, -1*(float)maxshift,"lag","l", n1, d1, o1,"sample index","i");
	    /* write to file */
    	sf_floatwrite(accumulate,n1*(2*maxshift+1),_accum);
    }	
	/* free accumulation arrays */
	free ( accumulate  );
	
	/* allocate array for warped trace */
	float* warped = sf_floatalloc(n1);
    /* apply shifts to the matching signal */
	dtw_apply_shifts( match, shifts, warped, n1);
	/* write the warped array to output */
    sf_floatwrite(warped,n1,_out);
	/* deallocate unneeded arrays */
	free (   match     );
	free (    ref      );
	free (   warped    );



	/* are we writing out the shifts? */	
	sf_file _shifts;
	float* shifts_out = NULL;
    if ( NULL != sf_getstring("shifts") ) {
		/* output integer shifts as floats */
		shifts_out = sf_floatalloc(n1);
		for ( i = 0 ; i < n1 ; i ++ ){
			shifts_out[i] = (float)shifts[i];
		}
	/* accumulation errors from forward and backtracking */ 
	    _shifts = sf_output ("shifts");
		/* set file parameters */
		dtw_set_file_params(_shifts, n1, d1, o1,"sample index","i", 1, 1, 0,"","");
	    /* write to file */
    	sf_floatwrite(shifts_out,n1,_shifts);
    }
	/* free shifting arrays */
	free ( shifts_out  );
	free (   shifts    );

    /* exit program */
    exit (0);
}
