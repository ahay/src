/* frequency response of PWD operator */

#include <rsf.h>
#include "lpwd.h"
#include "opwd.h"

int main(int argc, char*argv[])
{
	sf_file out;
	int n1, nf, n3, i3, interp;
	sf_complex **buf, c1;
	float d3, o3, radius;
	bool iir, opwd;

	sf_init(argc, argv);

	out = sf_output("out");

    if(!sf_getint("interp", &interp)) interp=0;
    /* 0: maxflat; 1: lagrange */
	if(!sf_getint("n1", &n1)) n1=50;
	/* samples in frequency domain between (0:f_c] */
	if(!sf_getint("nf", &nf)) nf=1;
	/* order of PWD */
	if(!sf_getfloat("o3", &o3)) o3=20;
	/* first dip angle */
	if(!sf_getfloat("d3", &d3)) d3=30;
	/* dip angle increment */
	if(!sf_getint("n3", &n3)) n3=1;
	/* number dip angle */

	sf_putint(out, "n1", 2*n1+1);
	sf_putfloat(out, "o1", -0.5);
	sf_putfloat(out, "d1", 0.5/n1);
	sf_putint(out, "n2", 2*n1+1);
	sf_putfloat(out, "o2", -0.5);
	sf_putfloat(out, "d2", 0.5/n1);
	sf_putint(out, "n3", n3);
	sf_putfloat(out, "d3", d3);
	sf_putfloat(out, "o3", o3);

	buf = sf_complexalloc2(2*n1+1, 2*n1+1);
	sf_settype(out,SF_COMPLEX);

	if(!sf_getbool("iir", &iir)) iir=false;
	/* y: iir; n: fir */

	if(!sf_getbool("opwd", &opwd)) opwd=true;
	/* y: circle interpolating pwd; n: line interpolating pwd */
	if(!sf_getfloat("radius", &radius)) radius=1.0;
	/* radius for circle interpolating pwd */

	if(opwd==true)
	{
		opwd_init(interp, nf);
		for(i3=0; i3<n3; i3++)
		{
			c1 = sf_cmplx(0, (d3*i3+o3)/180*SF_PI);
			opwd_freq(radius*cexpf(c1), n1, buf, iir);
			sf_complexwrite(buf[0], (2*n1+1)*(2*n1+1), out);
		}
		opwd_close();
	}else{
		lpwd_init(interp, nf, 0, 0);
		for(i3=0; i3<n3; i3++)
		{
			lpwd_freq((d3*i3+o3)/180*SF_PI, n1, buf, iir);
			sf_complexwrite(buf[0], (2*n1+1)*(2*n1+1), out);
		}
		lpwd_close();
	}
	free(buf[0]);
	free(buf);

	return 0;
}


