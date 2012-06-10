/* phase response of linear phase approximators */

#include <rsf.h>
#include "lphmf.h"
#include "lphlag.h"
#include "lphbsp.h"

int main(int argc, char*argv[])
{
	sf_file out;
	int n1, n2, nf, interp, i1, i2;
	sf_complex *obuf;
	float **bf, d2, o2, **c;

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


	sf_putint(out, "n1", n1);
	sf_putfloat(out, "o1", 0.0);
	sf_putfloat(out, "d1", 1.0/(n1-1));
	sf_putint(out, "n2", n2);
	sf_putfloat(out, "o2", o2);
	sf_putfloat(out, "d2", d2);

	bf = sf_floatalloc2(n1, n2);
	obuf = sf_complexalloc(2*n1-1);


	switch(interp)
	{
	case 1:
		c = lphlag(nf);
		for(i2=0; i2<n2; i2++)
		{
			lphlag_freq(nf, c, o2+d2*i2, n1-1, obuf);
			for(i1=0; i1<n1; i1++)
				bf[i2][i1] = cargf(obuf[i1+n1-1]);
		}
		break;
	case 2:
		c = lphbsp(nf);
		for(i2=0; i2<n2; i2++)
		{
			lphbsp_freq(nf, c, o2+d2*i2, n1-1, obuf);
			for(i1=0; i1<n1; i1++)
				bf[i2][i1] = cargf(obuf[i1+n1-1]);
		}
		break;
	default:
		c = lphmf(nf);
		for(i2=0; i2<n2; i2++)
		{
			lphmf_freq(nf, c, o2+d2*i2, n1-1, obuf);
			for(i1=0; i1<n1; i1++)
				bf[i2][i1] = cargf(obuf[i1+n1-1]);
		}
	}


	sf_floatwrite(bf[0], n1*n2, out);
	free(bf[0]);
	free(bf);
	free(c[0]);
	free(c);
	free(obuf);
	return 0;
}




