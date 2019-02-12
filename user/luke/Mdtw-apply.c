/* 
   program applies integer shifts (stored as floats!) to warp a matching trace.  Can match 1d shifts to a 1,2, or 3d volume , or 2d shifts to a 2, or 3 d volume, or 3d shifts to a 3d volume 
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
	int multishft = 0;
	int sn1, sn2, sn3;
    if (!sf_histint  (_shifts,"n1",&sn1))   sf_error("No n1= in shifts");
	if ( sn1 != n1 ) sf_error("shifts have different n1 than input");
    if (!sf_histint  (_shifts,"n2",&sn2))   sn2 = 1;
    if (!sf_histint  (_shifts,"n3",&sn3))   sn3 = 1;
	/* are we just using one shift ?*/
	if (sn3 == n3) multishft += 1;
	if (sn2 == n2) multishft += 2;
	
	int i2, i3;
    /* allocate input matching trace array */
    float* match  = sf_floatalloc(n1);
    /* allocate array for warped trace */
    float* warped = sf_floatalloc(n1);
    /* convert shifts to integer */
    int* shifts = sf_intalloc(n1);
    /* allocate secondary input reference shifts array*/
    float* shifts_f     = sf_floatalloc(n1);
	/* optional 2d shifts array */
	float* large_shifts ;
	/* and read from file if needed */
	if (multishft == 2){
		large_shifts = sf_floatalloc(n1*n2);
		/* case of 2d shifts for 3d volume */
		sf_floatread( large_shifts, n1*n2, _shifts);	
	}
	if (multishft == 0){
		/* case of only first axis matching */
		sf_floatread( shifts_f, n1, _shifts);
        /* convert float shifts to integer */
        for ( i = 0 ; i < n1 ; i++){
	        shifts[i] = (int)shifts_f[i];
        }
	}
	for (i3 = 0; i3 < n3 ; i3++){
		for ( i2 = 0; i2< n2 ; i2++){
			/* read shifts if applicable */
			if ( multishft >= 2 ){
			    if ( multishft == 3 ) {
				    /* everything matching */
				    sf_floatread( shifts_f, n1, _shifts);
			    } else {
			    	/* pull this shift from 2d shift array */
					shifts_f = _dtw_get_column( large_shifts, i2, n1 );
			    }
	            /* convert float shifts to integer */
	            for ( i = 0 ; i < n1 ; i++){
		            shifts[i] = (int)shifts_f[i];
	            }
		    }
	        /* and read from file */
            sf_floatread(match,n1,_in);
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
	if (multishft == 2){
		free (large_shifts);
	}

    /* exit program */
    exit (0);
}
