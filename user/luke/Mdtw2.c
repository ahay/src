/* 
   program warps a 2D input image to a 2D reference image
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
    if ( NULL == sf_getstring("ref") ) { sf_error("Need reference image ref=") ;}

    _ref  = sf_input ("ref");
	/* reference trace */ 
    _out = sf_output("out");
	int n1, n2, maxshift;
	float d1, d2, o1, o2, ex, str1, str2;
    /* Get sampling info */
    if (!sf_histint  (_in,"n1",&n1))   sf_error("No n1=");
    if (!sf_histfloat(_in,"d1",&d1))   sf_error("No d1=");
	if (!sf_histfloat(_in,"o1",&o1))   sf_error("No o1=");
    if (!sf_histint  (_in,"n2",&n2))   sf_error("No n2=");
    if (!sf_histfloat(_in,"d2",&d2))   sf_error("No d2=");
	if (!sf_histfloat(_in,"o2",&o2))   sf_error("No o2=");
    /* maximum shift */
    /* Get smoothing radius */
    if (!sf_getint("maxshift",&maxshift)){   
		sf_warning("maxshift set to 20");
		maxshift=20;
	}
    /* maximum shift */
    if (!sf_getfloat("exp",&ex))   ex = 2;
    /* error exponent (g-f)^exp */
	if (ex < 1){ 
		sf_warning("Exponent doesn't define a norm, changing to exp=2");
		ex = 2;
	}
    if (!sf_getfloat("strain1",&str1))   str1 = 1.0;
    /* maximum strain in first axis */
	if (str1 > 1){ 
		sf_warning("Program not set up for du/dt > 1, changing to str = 1");
		str1 = 1.0;
	}
    if (!sf_getfloat("strain2",&str2))   str2 = 1.0;
    /* maximum strain in second axis, if greatrer than 1 no strain limit */

	/* build N array for transpose operations */
	int ndim = 3;
	int* N = sf_intalloc(ndim);
	N[0] = 2*maxshift+1;
	N[1] = n1;
	N[2] = n2;
	
	/* build strain array */
	float* S = sf_floatalloc(ndim-1);
	S[0] = str1;
	S[1] = str2;
	
	/* number of alternations for smoothin */
	int nalter = 2;
	
	int fullsize = (int)dtw_size(N,ndim);
	int imgsize  = fullsize/N[0];
    /* allocate input matching image array */
    float* match  = sf_floatalloc(imgsize);
	/* and read from file */
	sf_floatread(match,imgsize,_in);
	/* allocate secondary input reference trace array*/
	float* ref    = sf_floatalloc(imgsize);
	/* and read from file */
    sf_floatread( ref,imgsize,_ref);
	/* decleare mismatch array */
	float* mismatch = sf_floatalloc(fullsize);
	/* declare traces for matching,  reference, and mismatch */
	float* matchtr    ;
	float* reftr      ;
	float* mismatchtr ;
    /* loop through traces and calculate mismatch */
    int i2;
	for ( i2 = 0 ; i2 < imgsize/N[1] ; i2++ ){
		/* get traces */
		matchtr = _dtw_get_column( match, i2, N[1] );
		reftr   = _dtw_get_column( ref  , i2, N[1] );
		/* calculate mismatch */		
		mismatchtr = _dtw_alignment_errors( matchtr, reftr, N[1], maxshift, ex);
		/* spread over null values */
		dtw_spread_over_nulls( mismatchtr, N[1], maxshift);
		/* write to global mismatch array */
        dtw_put_column( mismatch, mismatchtr, i2, N[0]*N[1] );
	}
	free (mismatchtr) ;
	free (reftr);
	free (matchtr);
    /* perform accumulation */
	int ismth;

    /* backup mismatch for backtracking */
	float* bu_mismatch = sf_floatalloc(fullsize);
	
    for ( ismth = 0; ismth < 2*nalter ; ismth ++){
    	/* loop through midpoints */
		    for ( i2 = 0 ; i2 < imgsize/N[1] ; i2++ ){
		    /* get mismatch */
			mismatchtr = _dtw_get_column( mismatch, i2, N[0]*N[1]) ;
			/* accumulate error */
	        mismatchtr = dtw_symmetric_accumulation(mismatchtr, N[1], maxshift, S[0]);
			/* subtract minimum error accumulation value */
			mismatchtr = dtw_subtract_minimum( mismatchtr, N[1]*N[0]);
			/* put back in array */
	    	dtw_put_column( mismatch, mismatchtr, i2, N[0]*N[1] );
	   	}
		/* determine if we are doing multiple smoothings */
		if (S[1] > 1) break;
		/* do we need to backup ?*/
		if ( ismth == 2*nalter - 2){
			dtw_acopy(bu_mismatch, mismatch, fullsize);
		}
	   	/* transpose */
    	mismatch = dtw_transp( mismatch, 2, 3, N, ndim);
	   	/* switch N */
	   	N = dtw_swaperoo(N, 2, 3, ndim);
		/* switch S */
		S = dtw_swaperoof(S, 1, 2, ndim-1);
	   	free (mismatchtr);
	}
	
	/* are we writing out the accumulation errors? */
	sf_file _accum;
    if ( NULL != sf_getstring("accum") ) {
		_accum = sf_output("accum");
		/* set parameters */
    	sf_putint   (_accum,"n1",N[0]); 
	    sf_putfloat (_accum,"d1",1);
	    sf_putfloat (_accum,"o1",-1*(float)maxshift);	
	    sf_putint   (_accum,"n2",N[1]); 
	    sf_putfloat (_accum,"d2",d1);
	    sf_putfloat (_accum,"o2",o1);
	    sf_putint   (_accum,"n3",N[2]);
	    sf_putfloat (_accum,"d3",d2);
	    sf_putfloat (_accum,"o3",o2);
		/* write accumulation errors */
		sf_floatwrite(mismatch,fullsize,_accum);
	}
	
	/* now lets figure out our image shifts */

	/* declare accumulation trace array */
	float* acumtr;
	/* declare shifts array */
	int* shifts = sf_intalloc(N[1]);
	/* allocate warped array */
	float* warped = sf_floatalloc(N[1]);
	
	/* set up optional shifts output */
	sf_file _shifts;
	/* and the switch */
	int shswitch ;
    if ( NULL != sf_getstring("shifts") ) {
		_shifts = sf_output("shifts");
		/* set parameters */
    	sf_putint   (_shifts,"n1",N[1]); 
	    sf_putfloat (_shifts,"d1",d1);
	    sf_putfloat (_shifts,"o1",o1);	
	    sf_putint   (_shifts,"n2",N[2]); 
	    sf_putfloat (_shifts,"d2",d2);
	    sf_putfloat (_shifts,"o2",o2);
		shswitch = 1;
	}	
	for ( i2 = 0 ; i2 < imgsize/N[1] ; i2++ ){
	    /* get mismatch */
		acumtr     = _dtw_get_column(    mismatch, i2, N[0]*N[1]) ;
		mismatchtr = _dtw_get_column( bu_mismatch, i2, N[0]*N[1]) ;
		/* backtrack */
		dtw_backtrack( acumtr, mismatchtr, shifts, N[1], maxshift, S[0]);
		/* get match trace */
		matchtr = _dtw_get_column( match, i2, N[1] );
	    /* apply shifts to the matching signal */
		dtw_apply_shifts( matchtr, shifts, warped, N[1]);
		/* write the warped signal to file */
	    sf_floatwrite(warped,N[1],_out);
		/* writing shifts?*/
		if (shswitch == 1){
		    sf_floatwrite(dtw_int_to_float( shifts, N[1]),N[1],_shifts);
	    }
	}
	
    /* free up arrays */	
    free ( mismatch );
	free ( bu_mismatch);
	free ( shifts );
	free ( warped );
	free ( matchtr );
	free ( acumtr );
	free ( ref );
	free ( match );

    /* exit program */
    exit (0);
}
