/* phase response of linear phase approximators */

#include <rsf.h>
#include "lphase.h"
#include "dsp.h"

int main(int argc, char*argv[])
{
	sf_file out;
	int n1, n2, nf, interp, i1, i2;
	sf_complex z1;
	float **bf, o1, d1, o2, d2, **c, *b1, f;

	sf_init(argc, argv);

	out = sf_output("out");

	if(!sf_getint("interp", &interp)) interp=0;
	/* 0: maxflat; 1: lagrange; 2: b-spline */
	if(!sf_getint("n1", &n1)) n1=50;
	/* samples in frequency domain between (0:f_c] */
	if(!sf_getint("nf", &nf)) nf=1;
	/* order of filters */
	if(!sf_getfloat("o2", &o2)) o2=0.1;
	/* first phase shift */
	if(!sf_getfloat("d2", &d2)) d2=0.3;
	/* phase shift increment */
	if(!sf_getint("n2", &n2)) n2=1;
	/* number of phase shift */

	o1 = 0.0;
	d1 = 0.5/(n1-1);
	sf_putint(out, "n1", n1);
	sf_putfloat(out, "o1", o1);
	sf_putfloat(out, "d1", d1);
	sf_putint(out, "n2", n2);
	sf_putfloat(out, "o2", o2);
	sf_putfloat(out, "d2", d2);

	bf = sf_floatalloc2(n1, n2);
	b1 = sf_floatalloc(2*nf+1);

	c = lphase(nf, interp);
	for(i2=0; i2<n2; i2++)
	{
		lphase_filt(nf, c, o2+d2*i2, b1, 2*nf, false);
		for(i1=0; i1<n1; i1++)
		{
			f = (o1+d1*i1);
			switch(interp)
			{
			case 1:
				z1 = fir_freq(-nf, nf, b1+nf, f);
				break;
			case 2:
				z1 = iir_freq(-nf, nf, b1+nf, -nf, nf, c[0]+nf, f);
				break;
			default:
				z1 = allpass_freq(nf, b1+nf, f);
			}
			bf[i2][i1] = cargf(z1);
		}
	}

	sf_floatwrite(bf[0], n1*n2, out);
	free(bf[0]);
	free(bf);
	free(c[0]);
	free(c);
	free(b1);
	return 0;
}




