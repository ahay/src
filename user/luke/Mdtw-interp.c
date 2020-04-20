/* 
   program takes traces sampled at arbitrary locations along a 1d line 
   and outputs regularly sampled line.
   program assumes the locations are ordered.
*/

#include <rsf.h>
#include "dtw.h"
int main (int argc, char* argv[])
{
	/* declare files */
    sf_file _in, _loc, _out;
	/* initialize rsf */
    sf_init (argc,argv);
	/* get files from inputs */
    _in  = sf_input("in");
	/* locations along line for initial sampling */
    if ( NULL == sf_getstring("loc") ) { sf_error("Need trace locations along line loc=") ;}
	_loc = sf_input ("loc");
	/* input sampling */
	int in_n1, in_n2, maxshift;
	float in_d1, in_o1;
	/* input location info */
	int loc_n1; 
    /* get sampling info from traces */
    if (!sf_histint  (_in,"n1",&in_n1))   sf_error("No n1= in input traces");
    if (!sf_histfloat(_in,"d1",&in_d1))   sf_error("No d1= in input traces");
	if (!sf_histfloat(_in,"o1",&in_o1))   sf_error("No o1= in input traces");
    if (!sf_histint  (_in,"n2",&in_n2))   sf_error("No n2= in input traces");
	/* does the location file have locations for each trace? */
	if (!sf_histint  (_loc,"n1",&loc_n1))   sf_error("No n1= in input trace locations");
	if (loc_n1 != in_n2 ) sf_error("Number of traces doesn't match number of locations.\nNeed n2 of input file = n1 of loc file.");
	/* dtw parameters*/
	float ex, str;
    /* maximum shift */
    if (!sf_getint("maxshift",&maxshift)){ 
		/* maximum shift to be tested */  
		sf_warning("maxshift set to 20");
		maxshift=20;
	}
    if (!sf_getfloat("exp",&ex))   ex = 2;
    /* error exponent (g-f)^exp */
	if (ex < 1){ 
		sf_warning("Exponent doesn't define a norm, changing to exp=2");
		ex = 2;
	}
    if (!sf_getfloat("strain",&str))   str = 1.0;
    /* maximum strain for dtw*/
	if (str > 1){ 
		sf_warning("Program not set up for du/dt > 1, changing to str = 1");
		str = 1.0;
	}
	/* allocate the input locations file */
	float* loc = malloc(sizeof(float)*loc_n1);
	/* read input locations file */
	sf_floatread( loc,loc_n1, _loc);
	/* get output sampling info	*/
	int out_n1, out_n2;
	float out_d1, out_d2, out_o1, out_o2 ;
	/* the output n1, o1, d1 are equal to the input since they are the traces */
	out_n1 = in_n1;
	out_d1 = in_d1;
	out_o1 = in_o1;	
	/* if nothing given, we will regularly resample */
	if (!sf_getint("n",&out_n2)){
		/* number of traces in output */
		out_n2 = in_n2;
		if (out_n2 < 2){
			sf_warning("Need more than 1 output trace, setting n=n2 of input");
			out_n2 = in_n2;
		}
	}
	/* lets make it so we are doing interior interpolation, keeping the first and last traces */
	out_o2 = loc[0];
	/* and lets determine our sampling */
	out_d2 = (loc[loc_n1-1] - loc[0])/(out_n2-1);
	/* let's make an array of output locations */
	float* places = malloc(sizeof(float)*out_n2);
	/* let's assign the place values */
	for (int i = 0 ; i < out_n2 ; i++){
		places[i] = out_o2 + i*out_d2;
	}
	/* prepare the output file */
    _out = sf_output("out");	
	/* write its sampling info to header file */
	sf_putint   (_out,"n1",out_n1); 
    sf_putfloat (_out,"d1",out_d1);
    sf_putfloat (_out,"o1",out_o1);	
    sf_putint   (_out,"n2",out_n2); 
    sf_putfloat (_out,"d2",out_d2);
    sf_putfloat (_out,"o2",out_o2);
		/* lets get the prior and next traces intialized */
	float* prior = malloc(sizeof(float)*in_n1);
	float* next  = malloc(sizeof(float)*in_n1);
	/* and read them from disk */
	sf_floatread( prior, in_n1, _in);
	sf_floatread( next , in_n1, _in);
	/* get their corresponding locations */
	float loc_prior = loc[0];
	float loc_next  = loc[1];
	/* and an integer to track where the next array is */
	int loc_pos = 1;
	/* let's initalize the interpolated trace */
	float* interpolated = malloc(sizeof(float)*out_n1);
	/* arrays for holding dtw shifts */
	int* shifts_forward = malloc(sizeof(int)*in_n1);
	int* shifts_back    = malloc(sizeof(int)*in_n1);
	float* shifts_fwd_float = malloc(sizeof(float)*in_n1);
	float* shifts_bk_float = malloc(sizeof(float)*in_n1);
	/* shifted traces for interpolation */
	float* shifted_prior = malloc(sizeof(int)*in_n1);
	float* shifted_next  = malloc(sizeof(int)*in_n1);
	/* initially say our shifts are not updated */
	bool updated_shifts = false;
	/* now lets start looping through our output locations */
	for ( int i = 0 ; i < out_n2 ; i++ ){
		/* where is our interpolated trace ?*/
		float spot = places[i];
		/* are we reading new traces this iteration? */
		/* check to make sure we are "sandwitched" by the right traces*/
		while ( !(loc_prior <= spot && loc_next >= spot) || (loc_pos >= in_n2 - 1) ){
			/* need new traces */
			updated_shifts = false;
			/* copy the next to prior */
			dtw_acopy(prior, next, in_n1);
			/* copy positions across */
			loc_prior = loc_next ;
			/* increment position */
			loc_pos++;
			/* get the next trace */
			sf_floatread (next, in_n1, _in);
			/* get the next location position */
			loc_next = loc[loc_pos];
		}
		/* are we on top of one of the traces ?*/
		if ( loc_prior == spot ){
			/* copy prior array to interpolated */
			dtw_acopy(interpolated, prior, out_n1);	
		}else{
			/* check if the next trace is in the interpolated spot */
			if ( loc_next == spot ){
				/* copy next array to interpolated */
				dtw_acopy(interpolated, next, out_n1);
			} else {
				/* do the dtw interpolation */
				/* how far apart are the ref and matching traces ? */
				float cum_dst = loc_next - loc_prior;
				/* how far forward are we ?*/
				float fwd_dst = spot - loc_prior ;
				/* and what is the relative fraction forward ? */
				float fwd_frac = fwd_dst / cum_dst ; 
				/* do we need to update our shifts ? */
				if (!updated_shifts){
					/* find shifts mapping the prior trace to the next if new traces */
					dtw_find_shifts( shifts_forward, next, prior, in_n1, maxshift, str, ex);
					/* find shifts mapping the next trace to the prior if new traces */
					dtw_find_shifts( shifts_back, prior, next, in_n1, maxshift, str, ex);
					/* and now we've updated our shifts! */
					updated_shifts = true;
				}
				/* warp the prior trace to our position */
				dtw_apply_scaled_shifts( prior, shifts_forward, shifted_prior, fwd_frac, in_n1);
				/* warp the next trace to our position */
				dtw_apply_scaled_shifts( next, shifts_back, shifted_next,1.0-fwd_frac, in_n1);	
				/* combine to create the interpolated trace */
				dtw_linear_comb( interpolated, 1.0-fwd_frac, shifted_prior, fwd_frac, shifted_next, in_n1);	
			}
		}
		/* write the interpolated trace */
		sf_floatwrite(interpolated, out_n1, _out);
	}
	/* free locations array */
	free ( loc) ;
	/* free output places array */
	free ( places );
	/* free traces used in interpolation */
	free ( prior );
	free ( next );
	free ( interpolated );
	/* free mapping shifts */
	free ( shifts_forward );
	free ( shifts_back );
	free ( shifts_fwd_float );
	free ( shifts_bk_float );
	/* free shifted traces */
	free ( shifted_prior );
	free ( shifted_next );
	/* remove temporary files */
	sf_close();
	/* program's over */
	exit (0);
}
