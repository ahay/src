/* 2-D smoothing by Gaussian filtering.

Takes: < input.rsf > smooth.rsf
*/

#include <rsf.h>

#include "gauss2.h"
#include "freqfilt2.h"

int main(int argc, char* argv[])
{
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

    if (!sf_getfloat("rect1",&rect1)) rect1=1.;
    if (!sf_getfloat("rect2",&rect2)) rect2=1.;
    /* smoothing radius */

    gauss2_init(n1, n2, rect1, rect2);

    input = sf_floatalloc(n12);
    smooth = sf_floatalloc(n12);

    for (i3=0; i3 < n3; i3++) {
	sf_read(input,sizeof(float),n12,in);
	
	freqfilt2_lop(false,false,n12,n12,input,smooth);

	sf_write(smooth,sizeof(float),n12,out);
    }

    exit(0);
}
