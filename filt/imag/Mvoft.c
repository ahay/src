/* V(t) function for a linear V(Z) profile.

Takes: < in.rsf > voft.rsf
*/

#include <math.h>

#include <rsf.h>

int main(int argc, char** argv)
{
    int n1, n2, i1, i2;
    float* vint;
    float o1, d1, v0, alpha, t;
    sf_file in, out;

    sf_init (argc, argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint (in,"n1",&n1)) sf_error ("No n1= in input\n");
    if (!sf_histint (in,"n2",&n2)) sf_error ("No n2= in input\n");

    if (!sf_histfloat (in,"o1",&o1)) o1=0.;
    if (!sf_histfloat (in,"d1",&d1)) sf_error ("No d1= in input\n");

    if (!sf_getfloat ("v0",&v0)) v0 = 1.5;
    /* initial velocity */
    if (!sf_getfloat ("alpha",&alpha)) alpha = 0.5;
    /* velocity gradient */

    vint = sf_floatalloc(n1);
    for (i1=0; i1 < n1; i1++) {
	t = alpha*(o1+i1*d1);
	if (t > 0.) {
	    vint[i1] = v0 * sqrtf((expf(t)-1.)/t);
	} else {
	    vint[i1] = v0;
	}
    }
    for (i2=0; i2 < n2; i2++) {
	sf_write(vint,sizeof(float),n1,out);
    }

    sf_close();
    exit (0);
}

/* 	$Id: Mvoft.c,v 1.3 2004/03/22 05:43:24 fomels Exp $	 */

