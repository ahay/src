/* frequency response of PWD operator */

#include <rsf.h>
#include "pwd.h"

int main(int argc, char*argv[])
{
	sf_file out;
	int n1, nw, n3, i;
	sf_complex **d;
	float d3, o3, r;
	bool iir, opwd;

	sf_init(argc, argv);

	out = sf_output("out");

	if(!sf_getint("n1", &n1)) n1=50;
	/* samples in frequency domain between (0:f_c] */
	if(!sf_getint("nw", &nw)) nw=1;
	/* order of PWD */
	if(!sf_getfloat("o3", &o3)) o3=10;
	/* first dip angle */
	if(!sf_getfloat("d3", &d3)) d3=20;
	/* dip angle increment */
	if(!sf_getint("n3", &n3)) n3=1;
	/* number dip angle */
	if(!sf_getbool("iir", &iir)) iir=false;
	/* y: iir; n: fir */
	if(!sf_getbool("opwd", &opwd)) opwd=true;
	/* y: circle interpolating pwd; n: line interpolating pwd */
	if(!sf_getfloat("r", &r)) r=1.0;
	/* radius for opwd */
	
	d = sf_complexalloc2(2*n1+1, 2*n1+1);

	sf_settype(out,SF_COMPLEX);
	sf_putint(out, "n1", 2*n1+1);
	sf_putfloat(out, "o1", -1.0);
	sf_putfloat(out, "d1", 1.0/n1);
	sf_putint(out, "n2", 2*n1+1);
	sf_putfloat(out, "o2", -1.0);
	sf_putfloat(out, "d2", 1.0/n1);
	sf_putint(out, "n3", n3);
	sf_putfloat(out, "d3", d3);
	sf_putfloat(out, "o3", o3);

	pwd_init(nw,0,0);
	for(i=0; i<n3; i++)
	{
		if(opwd) opwd_freq(	r*sin((d3*i+o3)/180*M_PI), 
					r*cos((d3*i+o3)/180*M_PI), n1, d, iir);
		else 	lpwd_freq(tan((d3*i+o3)/180*M_PI), n1, d, iir);
		sf_complexwrite(d[0], (2*n1+1)*(2*n1+1), out);
	}

	pwd_release();
	return 0;
}


