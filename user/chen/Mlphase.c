/* phase response of linear phase approximators */

#include <rsf.h>
#include "lphase.h"

int main(int argc, char*argv[])
{
	sf_file out;
	int n1, n2, nw, interp, i1, i2;
	sf_complex *ph;
	float **bf, d2, o2;
	void *h;
	bool causal;

	sf_init(argc, argv);

	out = sf_output("out");

	if(!sf_getint("interp", &interp)) interp=0;
	/* 0: maxflat; 1: lagrange */
	if(!sf_getint("n1", &n1)) n1=50;
	/* samples in frequency domain between (0:f_c] */
	if(!sf_getint("nw", &nw)) nw=1;
	/* order of filters */
	if(!sf_getfloat("o2", &o2)) o2=10;
	/* first dip angle */
	if(!sf_getfloat("d2", &d2)) d2=20;
	/* dip angle increment */
	if(!sf_getint("n2", &n2)) n2=1;
	/* number dip angle */
	if(!sf_getbool("causal", &causal)) causal=false;
	/* y: causal; n: noncausal */

	h = lphase_init(nw, causal, interp);

	sf_putint(out, "n1", n1);
	sf_putfloat(out, "o1", 0.0);
	sf_putfloat(out, "d1", 1.0/n1);
	sf_putint(out, "n2", n2);
	sf_putfloat(out, "o2", o2);
	sf_putfloat(out, "d2", d2);

	bf = sf_floatalloc2(n1, n2);
	ph = sf_complexalloc(2*n1+1);

	for(i2=0; i2<n2; i2++)
	{
		lphase_freq(h, tan((d2*i2+o2)/180*SF_PI), n1, ph);
		for(i1=0; i1<n1; i1++)
			bf[i2][i1] = cargf(ph[i1+n1]);
	}

	sf_floatwrite(bf[0], n1*n2, out);
	lphase_close(h);
	free(ph);
	free(bf[0]);
	free(bf);
	return 0;
}


