/* Pad and interleave traces. 

Takes < small.rsf mask=mask.rsf > large.rsf

Each initial trace is followed by I<jump> zero traces, the same for planes.
*/

#include <rsf.h>

int main(int argc, char* argv[])
{
    int n1,n2,n3, jump, i1,i2,i3, j;
    float d, *pp, *zero, *one;
    sf_file in, out, mask;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    mask = sf_output("mask");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_getint("jump",&jump)) jump=2;
    /* aliasing */

    if (n2 > 1) {
	sf_putint(out,"n2",n2*jump);
	if (sf_histfloat(in,"d2",&d)) sf_putfloat(out,"d2",d/jump);
    }
    if (n3 > 1) {
	sf_putint(out,"n3",n3*jump);
	if (sf_histfloat(in,"d3",&d)) sf_putfloat(out,"d3",d/jump);
    }

    sf_fileflush(mask,out);

    pp = sf_floatalloc(n1);
    zero = sf_floatalloc(n1);
    one = sf_floatalloc(n1);

    for (i1=0; i1 < n1; i1++) {
	zero[i1] = 0.;
	one[i1] = 1.;
    }

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    sf_floatread(pp,n1,in);
	    sf_floatwrite(pp,n1,out);
	    sf_floatwrite(one,n1,mask);
	    if (n2 >1) {
		for (j=0; j < jump-1; j++) {
		    sf_floatwrite(zero,n1,out);
		    sf_floatwrite(zero,n1,mask);
		}
	    }
	}
	if (n3 > 1) {
	    for (i2=0; i2 < n2*jump; i2++) {
		sf_floatwrite(zero,n1,out);
		sf_floatwrite(zero,n1,mask);
	    }
	}
    }

    sf_close();
    exit(0);
}

