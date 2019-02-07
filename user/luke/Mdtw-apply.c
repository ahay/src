/* 
   program applies integer shifts (stored as floats!) to warp a matching trace 
*/

#include <rsf.h>
#include "dtw.h"

int main (int argc, char* argv[])
{
	/* declare files */
    sf_file _in, _shifts, _out;
	/* initialize rsf */
    sf_init (argc,argv);
	/* get files from inputs */
    _in = sf_input("in");
    if ( NULL == sf_getstring("shifts") ) { sf_error("Need shifts file shifts=") ;}

    _shifts  = sf_input ("shifts");
	/* reference trace */ 
    _out = sf_output("out");
	int n1, i;
	float d1, o1;
    /* Get sampling info */
    if (!sf_histint  (_in,"n1",&n1))   sf_error("No n1=");
    if (!sf_histfloat(_in,"d1",&d1))   sf_error("No d1=");
	if (!sf_histfloat(_in,"o1",&o1))   sf_error("No o1=");
    /* allocate input matching trace array */
    float* match  = sf_floatalloc(n1);
	/* and read from file */
	sf_floatread(match,n1,_in);
	/* allocate secondary input reference shifts array*/
	float* shifts_f    = sf_floatalloc(n1);
	/* and read from file */
    sf_floatread( shifts_f, n1, _shifts);
	/* convert shifts to integer */
	int* shifts = sf_intalloc(n1);
	/* convert float shifts to integer */
	for ( i = 0 ; i < n1 ; i++){
		shifts[i] = (int)shifts_f[i];
	}
	/* allocate array for warped trace */
	float* warped = sf_floatalloc(n1);
    /* apply shifts to the matching signal */
	dtw_apply_shifts( match, shifts, warped, n1);
	/* write the warped array to output */
    sf_floatwrite(warped,n1,_out);
	/* deallocate  arrays */
	free (   match     );
	free (   warped    );
	free (   shifts_f  );
	free (   shifts    );

    /* exit program */
    exit (0);
}
