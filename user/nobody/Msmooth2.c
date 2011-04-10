/* 2-D smoothing by Gaussian filtering.

Takes: < input.rsf > smooth.rsf
*/

#include <rsf.h>

#include "triangle2.h"
#include "gauss2.h"
#include "freqfilt2.h"

int main(int argc, char* argv[])
{
    bool gauss;
    int n1, n2, n12, n3, i3;
    float rect1, rect2, *input, *smooth;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    n12 = n1*n2;

    if (!sf_getfloat("rect1",&rect1)) rect1=3.;
    if (!sf_getfloat("rect2",&rect2)) rect2=3.;
    /* smoothing radius */

    if (!sf_getbool("gauss",&gauss)) gauss=false;
    /* if y, use exact Gaussian */
    
    if (gauss) {
	gauss2_init(n1, n2, rect1, rect2);
    } else {
	triangle2_init((int) rect1, (int) rect2, n1, n2);
    }

    input = sf_floatalloc(n12);
    smooth = sf_floatalloc(n12);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(input,n12,in);
	
	if (gauss) {
	    freqfilt2_lop(false,false,n12,n12,input,smooth);
	} else {
	    triangle2_lop(false,false,n12,n12,input,smooth);
	}

	sf_floatwrite(smooth,n12,out);
    }


    exit(0);
}

/* 	$Id$	 */

