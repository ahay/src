/* 
   problem finds the optimal trajectory 
   across accumulation errors using backtracking.  
   takes the same strain strain= as the accumulation step.
*/

#include <rsf.h>
#include "dtw.h"

int main (int argc, char* argv[])
{
	/* declare files */
    sf_file _in, _out;
	/* initialize rsf */
    sf_init (argc,argv);
	/* get files from inputs */
    _in = sf_input("in");
    _out = sf_output("out");
	
	int n1, n2, maxshift;
	float d1, d2, o1, o2, str;
    /* Get sampling info */
    if (!sf_histint  (_in,"n1",&n1))   sf_error("No n1=");
    if (!sf_histfloat(_in,"d1",&d1))   sf_error("No d1=");
	if (!sf_histfloat(_in,"o1",&o1))   sf_error("No o1=");
    if (!sf_histint  (_in,"n2",&n2))   sf_error("No n2=");
    if (!sf_histfloat(_in,"d2",&d2))   sf_error("No d2=");
	if (!sf_histfloat(_in,"o2",&o2))   sf_error("No o2=");
	
	/* convert n1 to maxshift */
	maxshift = (n1-1)/2;

    if (!sf_getfloat("strain",&str))   str = 1;
    /* maximum strain */
	if (str > 1){ 
		sf_warning("Program not set up for du/dt > 1, changing to str = 1");
		str = 1;
	}
	/* allocate the accumulation array */
	float* accumulate = sf_floatalloc(n1*n2);
	/* read the accumulation array from input */
	sf_floatread(accumulate,n1*n2,_in);
	/* allocate shifts array */
	int* shifts = sf_intalloc(n1);
	/* mismatch file used in backtracking if strain less than 1 */
	sf_file _miss;
	float* mismatch = sf_floatalloc(n1*n2);
	if (str < 1){
	    if ( NULL == sf_getstring("error") ) { sf_error("Need error file error= for strain < 1") ;}
	    _miss  = sf_input ("error");
		/* error or misfit to accumulate */
		sf_floatread(mismatch,n1*n2,_miss);
	}
	/* backtrack to find integer shifts */
	dtw_backtrack( accumulate, mismatch, shifts, n2, maxshift, str);
	/* convert the integer shifts to float for output */
	float* shifts_out = sf_floatalloc(n2);
    for (int i = 0 ; i < n2 ; i ++){
		shifts_out[i] = (float)shifts[i];
		}
	/* define output file parameters */ 
	sf_putint   (_out,"n1",n2);
	sf_putfloat (_out,"d1",d2);
	sf_putfloat (_out,"o1",o2);
	sf_putstring(_out,"label1", "index");
	sf_putstring(_out, "unit1",  "i");
	sf_putint   (_out,"n2",1);
	sf_putfloat (_out,"d2",1);
	sf_putfloat (_out,"o2",0);
	sf_putstring(_out,"label2", "");
	sf_putstring(_out, "unit2",  "");
	/* write to file */
    sf_floatwrite(shifts_out,n2,_out);
	/* free shifting arrays */
	free ( shifts_out  );
	free (   shifts    );
    free ( mismatch    );
	free ( accumulate  );
    /* exit program */
    exit (0);
}
