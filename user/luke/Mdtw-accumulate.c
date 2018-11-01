/* 
   accumulates or smooths alignment errors in the 
   forward direction, subject to 
   strain str=float is the maximum du/di
*/

#include <rsf.h>
#include "dtw.h"

int main (int argc, char* argv[])
{
	/* declare files */
    sf_file _in, _out;
	/* initialize rsf */
    sf_init (argc,argv);
	/* get alignment errors as input */
    _in = sf_input("in");
	/* prepare output file */
    _out = sf_output("out");
	int n1, n2, maxshift;
	float d1, d2, o1, o2, str;
    /* Get sampling info */
    if (!sf_histint  (_in,"n1",&n1))   sf_error("No n1=");
    if (!sf_histfloat(_in,"d1",&d1))   sf_error("No d1=");
	if (!sf_histfloat(_in,"o1",&o1))   sf_error("No o1=");
    if (!sf_histint  (_in,"n2",&n2))   sf_error("No n2=");   
	if (!sf_histfloat  (_in,"d2",&d2))   sf_error("No d2=");
	if (!sf_histfloat  (_in,"o2",&o2))   sf_error("No o2=");
	/* get strain info */
    if (!sf_getfloat("strain",&str))   str = 1;
    /* maximum strain */
	if (str > 1){ 
		sf_warning("Program not set up for du/dt > 1, changing to str = 1");
		str = 1;
	}
	int dir;
    if (!sf_getint("dir",&dir))   dir = 0;
    /* accumulation direction: 1 is forward, -1 is backward, 0 is both */
	
	/* determine max shift */
	maxshift = (int) (n1-1)/2 ;
    /* allocate input mismatch array */
	float* mismatch = sf_floatalloc(n1*n2);
	/* and read */
	sf_floatread(mismatch,n1*n2,_in);
	/* declare array for forward error accumulation */
	/* allocate the output array, accumulate */
	float* accumulate = sf_floatalloc(n1*n2);
	/* initialize the accumulation array, wit zeros if directional 
	   and with the negative errors if symmetric */
	if (dir == 0){
		dtw_acopy(accumulate,mismatch,n1*n2) ;
		dtw_mul(accumulate,-1., n1*n2) ;
	} else { 
		/* set to sero */
		dtw_copy(accumulate,0.,n1*n2);
	}
	float* accumulate_f ;
	/* see if we are accumulating forward */
		if (dir >= 0){ 
			/* allocate forward error array */
			accumulate_f = sf_floatalloc(n1*n2);
			/* accumulate forward errors */
			dtw_accumulate_errors( mismatch, accumulate_f, n2, maxshift, str);
			/* add the accumulation errors to the output array */
			dtw_aadd( accumulate, accumulate_f, n1*n2) ;
			/* free intermediate arrays */
			free ( accumulate_f );
		}

	/* declare array for backward error accumulation */
	float* accumulate_b ;
	/*  reverse mismatch */
	float* mismatch_r ;
	/* and the reversed reverse accumulation thats in the correct orientation */
	float* accumulate_br;
	/* now check to see if we are accumulating backward */
	if (dir <= 0){
		/* allocate arrays */
		accumulate_b  = sf_floatalloc(n1*n2);
		mismatch_r    = sf_floatalloc(n1*n2);
		accumulate_br = sf_floatalloc(n1*n2);
		/* reverse mismatch array */
		dtw_reverse_array(mismatch, mismatch_r, n2, n1);
		/* accumulate backward errors */
		dtw_accumulate_errors( mismatch_r, accumulate_b, n2, maxshift, str);
		/* blip the backward errors */
		dtw_reverse_array(accumulate_b,accumulate_br, n2, n1);
		/* add the backward accumulation to the output array */
		dtw_aadd( accumulate, accumulate_br, n1*n2) ;
		/* free intermediate arrays */
		free ( accumulate_b );
		free ( mismatch_r   );
		free ( accumulate_br);
	}

	/* write the accumulated error to output */
    sf_floatwrite(accumulate,n1*n2,_out);
	/* deallocate arrays */
	free (accumulate);
	free (mismatch)  ;


    exit (0);
}
