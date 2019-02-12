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
	
	/* see if the matching object is higher dimensional */
	int n2, n3;
	if (!sf_histint  (_in,"n2",&n2)) n2 = 1; 
	if (!sf_histint  (_in,"n3",&n3)) n3 = 1;
	int nmid = n2*n3;
	
	/* are the shifts also higher dimensional or compatable? */
	int multishft = 1;
	int sn1, sn2, sn3;
    if (!sf_histint  (_shifts,"n1",&sn1))   sf_error("No n1= in shifts");
	if ( sn1 != n1 ) sf_error("shifts have different n1 than input");
    if (!sf_histint  (_shifts,"n2",&sn2))   sn2 = 1;
    if (!sf_histint  (_shifts,"n3",&sn3))   sn3 = 1;
	/* are we just using one shift ?*/
	if (sn2*sn3 != nmid) multishft = 0;
	
	int i2, i3;
    /* allocate input matching trace array */
    float* match  = sf_floatalloc(n1);
    /* allocate array for warped trace */
    float* warped = sf_floatalloc(n1);
    /* convert shifts to integer */
    int* shifts = sf_intalloc(n1);
    /* allocate secondary input reference shifts array*/
    float* shifts_f    = sf_floatalloc(n1);
	
	/* and read from file */
	for (i3 = 0; i3 < n3 ; i3++){
		for ( i2 = 0; i2< n2 ; i2++){
			/* read shifts if applicable */
			if ( multishft == 1 || i2+i3 == 0) sf_floatread( shifts_f, n1, _shifts);
	        /* and read from file */
            sf_floatread(match,n1,_in);
	        /* convert float shifts to integer */
	        for ( i = 0 ; i < n1 ; i++){
		        shifts[i] = (int)shifts_f[i];
	        }
            /* apply shifts to the matching signal */
	        dtw_apply_shifts( match, shifts, warped, n1);
	        /* write the warped array to output */
            sf_floatwrite(warped,n1,_out);
		}
	}
	/* deallocate  arrays */
	free (   match     );
	free (   warped    );
	free (   shifts_f  );
	free (   shifts    );

    /* exit program */
    exit (0);
}
