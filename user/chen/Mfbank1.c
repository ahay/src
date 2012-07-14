/* 2d filter bank  */

#include <rsf.h>
#include "fbank.h"


int main(int argc, char*argv[])
{
	sf_file in, out;
	int nf, n1, n2, interp, nf2;
	int i2;
	float *wav, **fb;

	sf_init(argc, argv);

	in  = sf_input("in");
	out = sf_output("out");

	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

	if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
	n2 = sf_leftsize(in, 1);

	sf_shiftdim(in,out,1);

	if (!sf_getint("order",&nf)) nf=1;
	/* PWD filter order */
	
	if (!sf_getint("interp",&interp)) interp=0;
	/* interpolation method: 
	0: maxflat
	1: Lagrange 
	2: B-Spline */

	nf2 = 2*nf+1;
	wav = sf_floatalloc(n1);
	fb  = sf_floatalloc2(n1, nf2);

	sf_putint(out, "n2", nf2);
	sf_putfloat(out, "o2", 0);
	sf_putfloat(out, "d2", 1);


	fbank_init(nf, interp);


	for(i2=0; i2<n2; i2++)
	{
		sf_floatread(wav, n1, in);
		fbank1(n1, wav, fb);
		sf_floatwrite(fb[0], n1*nf2, out);
	}

	fbank_close();
	return 0;
}



