/* Apply 2-D directional high-pass to highlight data.

Takes: < inp.rsf > highlighted.rsf
*/

#include <rsf.h>

#include "hipass.h"

int main(int argc, char* argv[])
{
    int n1, n2, i1, i2;
    float ax, ay, eps, num, den;
    float **map, **dx, **dy, *t1, *t2;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");

    if (!sf_getfloat("ax",&ax)) ax=1.;
    /* x direction */
    if (!sf_getfloat("ay",&ay)) ay=1.;
    /* y direction */
    if (!sf_getfloat("eps",&eps)) eps=0.;
    /* highpass filter parameter; if eps=0, apply derivative */

    map = sf_floatalloc2(n1,n2);
    dx =  sf_floatalloc2(n1,n2);
    dy =  sf_floatalloc2(n1,n2);
    t1 = sf_floatalloc(n2);
    t2 = sf_floatalloc(n2);

    sf_floatread (map[0],n1*n2,in);

    if (eps > 0.) {
	hipass_init(eps);
	for (i2=0; i2 < n2; i2++) {
	    hipass_lop (false, false, n1, n1, map[i2], dx[i2]);
	}
	for (i1=0; i1 < n1; i1++) {
	    for (i2=0; i2 < n2; i2++) {
		t1[i2] = map[i2][i1];
	    }
	    hipass_lop (false, false, n2, n2, t1, t2);
	    for (i2=0; i2 < n2; i2++) {
		dy[i2][i1] = t2[i2];
	    }
	}
    } else {
	for (i2=0; i2 < n2; i2++) {
	    dx[i2][0] = 0.;
	    for (i1=1; i1 < n1-1; i1++) {
		dx[i2][i1] = map[i2][i1+1] - map[i2][i1-1];
	    }
	    dx[i2][n1-1] = 0.;
	}
	for (i1=0; i1 < n1; i1++) {
	    dy[0][i1] = 0.;
	    dy[n2-1][i1] = 0.;
	}
	for (i2=1; i2 < n2-1; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		dy[i2][i1] = map[i2+1][i1] - map[i2-1][i1];
	    }
	}
    }

    if (ax == 0. && ay==0.) {
	num = den = 0.;
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		num += dx[i2][i1]*dy[i2][i1];
		den += dy[i2][i1]*dy[i2][i1];
	    }
	}
	ax = 1.; 
	ay = - num/den;

	sf_warning("ay=%g",ay);
    }

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    map[i2][i1] = ax*dx[i2][i1]+ay*dy[i2][i1];
	}
    }

    sf_floatwrite (map[0],n1*n2,out);

    sf_close();
    exit(0);
}

