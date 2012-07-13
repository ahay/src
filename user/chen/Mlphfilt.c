/* linear phase filter */

#include <rsf.h>

#include "lphase.h"

int main(int argc, char* argv[])
{
	sf_file out;
	int nf, interp;
	float **c;

	sf_init(argc, argv);
	out = sf_output ("out");

	if(!sf_getint("order", &nf)) nf=1;
	/* order of linear phase filter */
	if (!sf_getint("interp",&interp)) interp=0;
	/* interpolation method: 
	0: maxflat
	1: Lagrange 
	2: B-Spline */

	sf_putint(out, "n1", 2*nf+1);
	sf_putfloat(out, "o1", -nf);
	sf_putfloat(out, "d1", 1);
	sf_putint(out, "n2", 2*nf+1);
	sf_putfloat(out, "o2", 0);
	sf_putfloat(out, "d2", 1);

	c = lphase(nf, interp);
	sf_floatwrite(c[0], (2*nf+1)*(2*nf+1), out);

	free(c[0]);
	free(c);
	return 0;
}

