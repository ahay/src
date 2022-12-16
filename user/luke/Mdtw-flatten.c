/* flattens a gather or similar object to its stack using dtw, optionally writes out shifts, currently set up for (time,gather,space) for 2d imaging
 */
#include <rsf.h>
#include "dtw.h"
int main (int argc, char* argv[])
{
    /* input gather files */
    sf_file _in, _shifts, _out;
    /* initialize rsf */
    sf_init (argc,argv);
    /* input gathers file */
    _in = sf_input("in");
    /* flattened gather file */ 
    _out = sf_output("out");


    /* declare sampling variables */
    int   n1, n2, n3;
    float d1, d2, d3;
    float o1, o2, o3;
	
    /* Get sampling info */
	
    /* time/depth sampling */
    if (!sf_histint  (_in,"n1",&n1))   sf_error("No n1=");
    if (!sf_histfloat(_in,"d1",&d1))   sf_error("No d1=");
    if (!sf_histfloat(_in,"o1",&o1))   sf_error("No o1=");
    /* gather sampling */
    if (!sf_histint  (_in,"n2",&n2))   sf_error("No n2=");
    if (!sf_histfloat(_in,"d2",&d2))   sf_error("No d2=");
    if (!sf_histfloat(_in,"o2",&o2))   sf_error("No o2=");	
    /* spatial sampling x */
    if (!sf_histint  (_in,"n3",&n3))   sf_error("No n3=");
    if (!sf_histfloat(_in,"d3",&d3))   sf_error("No d3=");
    if (!sf_histfloat(_in,"o3",&o3))   sf_error("No o3=");	


    /* shifts switch */
    int sh = 0; 
    if ( NULL != sf_getstring("shifts") ) {
	/* output gather flattening shifts */ 
	_shifts = sf_output ("shifts");
	sh = 1;
	/* format shifts file */
	sf_putint   (_shifts,"n1",n1); 
	sf_putfloat (_shifts,"d1",d1);
	sf_putfloat (_shifts,"o1",o1);
	sf_putint   (_shifts,"n2",n2);
	sf_putfloat (_shifts,"d2",d2);
	sf_putfloat (_shifts,"o2",o2);
	sf_putint   (_shifts,"n3",n3);
	sf_putfloat (_shifts,"d3",d3);
	sf_putfloat (_shifts,"o3",o3);
    } else {
	_shifts = NULL;
    }
	
    /* get dtw parameters */
	
    float ex;
    if (!sf_getfloat("exp",&ex))   ex = 2;
    /* error exponent (g-f)^exp */
    if (ex < 1){ 
	sf_warning("Exponent doesn't define a norm, changing to exp=2");
	ex = 2;
    }
    float str;
    if (!sf_getfloat("strain",&str))   str = 1.0;
    /* maximum strain */
    if (str > 1){ 
	sf_warning("Program not set up for du/dt > 1, changing to strain = 1");
	str = 1.0;
    }
    if ((int)((float)n1/10) < (int)(1/str) || str == 0){
	str = (1/((float)n1/10));
	sf_warning("Strain too small! Changing to strain = %g",str);
    }
    int maxshift;
    if (!sf_getint("maxshift",&maxshift)){ maxshift=20;
	/* maximum shift */
	sf_warning("maxshift set to 20");

    }
	
    int ig, ng, io, no;
    /* set number of gathers */
    ng = n3;
    /* number of offsets in gather */
    no = n2;
    /* declare null value */
    float nullval = -9999;
	
    /* declare gather input array */
    float* gather = malloc(sizeof(float)*n1*no);
    /* declare nulled gather array */
    float* ngather = malloc(sizeof(float)*n1*no);
    /* declare stack it will be matched to */
    float* stack  = malloc(sizeof(float)*n1);
    /* declare matching trace */
    float* match  = malloc(sizeof(float)*n1);
    /* declare match trace for nulls */
    float* nmatch  = malloc(sizeof(float)*n1);
    /* declare shifts array */
    int*  tr_shifts = malloc(sizeof(int)*n1);
    /* declare gather shifts */
    float* gather_shifts = malloc(sizeof(float)*n1*no);
    /* declare warped trace array */
    float* warped = malloc(sizeof(float)*n1);
    /* warped gather array */
    float* warped_gather = malloc(sizeof(float)*n1*no);

    /* loop through gathers g for gather, o for offsets*/
    for ( ig = 0; ig < ng ; ig++ ){
	/* read input gather */
	sf_floatread(gather,n1*no,_in);
	/* write null values */
	free (ngather);
	ngather = dtw_write_gather_nullvals(gather, n1, no, nullval);
	/* create stack for reference */
	dtw_norm_stack(ngather,stack,n1,no,nullval);
	/* and free */

	/* zero out warped gather */
	dtw_copy( warped_gather, 0., n1*no);
	/* zero out gather shifts */
	dtw_copy( gather_shifts, 0., n1*no);
	for ( io = 0 ; io < no ; io++){
	    /* get traces */
	    dtw_get_column( gather, match, io, n1);
	    /* determine shifts */
	    dtw_find_shifts( tr_shifts, stack, match, n1, maxshift, str, ex);
	    /* apply shifts to trace */
	    dtw_apply_shifts( match, tr_shifts, warped, n1);
	    /* write trace shifts to gather shifts array */
	    dtw_put_column( gather_shifts, dtw_int_to_float( tr_shifts, n1), io, n1 ) ; 
	    /* write warped trace to gather array */
	    dtw_put_column( warped_gather, warped, io, n1 )	; 
	}
	/* write out warped_gather */
	sf_floatwrite(warped_gather,n1*no,_out);
	/* write out shifts if desired */
	if (sh == 1){
	    sf_floatwrite(gather_shifts,n1*no,_shifts);
	}	
    }
	
    /* free arrays */
    free (  gather       );
    free (  ngather      );
    free (  stack        );
    free (  match        );
    free (  nmatch       );
    free (  tr_shifts    );
    free (  warped       );
    free ( gather_shifts );
    free ( warped_gather );

	
    /* exit program */
    exit (0);
}
