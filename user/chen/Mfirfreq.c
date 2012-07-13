/* frequency response of fir filters */

#include <rsf.h>

#include "dsp.h"

int main(int argc, char* argv[])
{
	sf_file in, out;
	int n1, n2, mw, o1, d1;
	int iw, i2;
	float **c, ow, dw;
	sf_complex **d;

	sf_init(argc, argv);
	in  = sf_input ("in");
	out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
	if (!sf_histint(in,"o1", &o1)) sf_error("o1");
	if (!sf_histint(in,"d1", &d1)) sf_error("d1");

	if(!sf_getint("nw", &mw)) mw=101;
	/* samples in frequency domain */
	if(!sf_getfloat("ow", &ow)) ow=-0.5;
	/* first frequency */
	if(!sf_getfloat("dw", &dw)) dw=0.01;
	/* frequency increment*/

	sf_putint(out, "n1", mw);
	sf_putfloat(out, "o1", ow);
	sf_putfloat(out, "d1", dw);
	sf_settype(out,SF_COMPLEX);

	c = sf_floatalloc2(n1, n2);
	d = sf_complexalloc2(mw, n2);

	sf_floatread(c[0], n1*n2, in);

	for(i2=0; i2<n2; i2++)
	for(iw=0; iw<mw; iw++)
	d[i2][iw] = fir_freq(o1, o1+d1*(n1-1), c[i2]-o1, (ow+dw*iw));

	sf_complexwrite(d[0], mw*n2, out);

	free(c[0]);
	free(c);
	free(d[0]);
	free(d);
	return 0;
}

