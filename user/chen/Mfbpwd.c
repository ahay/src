/* omnidirectional plane-wave destruction  */

#include <rsf.h>
#include "fbpwd.h"


int main(int argc, char*argv[])
{
	sf_file in, out, fp;
	int n1, n2, n3, m1, m2;
	int i3;
	float ****wav, **dip, radius;

	sf_init(argc, argv);

	in  = sf_input("in");
	out = sf_output("out");
	fp  = sf_input("dip");

	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");
	if (SF_FLOAT != sf_gettype(fp)) sf_error("Need float type");

	if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
	if (!sf_histint(in, "n2", &n2)) sf_error("No n2= in input");
	if (!sf_histint(in, "n3", &m1)) sf_error("No m1= in input");
	if (!sf_histint(in, "n4", &m2)) sf_error("No m2= in input");
	n3 = sf_leftsize(in, 4);
	sf_unshiftdim2(in,out,3);

	if (!sf_getfloat("radius", &radius)) radius = 1.0;
	/* interpolating radius for opwd */

	wav = sf_floatalloc4(n1, n2, m1, m2);
	dip = sf_floatalloc2(n1,n2);

	fbpwd_init(n1, n2, radius);

	for(i3=0; i3<n3; i3++)
	{
		sf_floatread(wav[0][0][0], n1*n2*m1*m2, in);
		sf_floatread(dip[0], n1*n2, fp);
		fbpwd(m1, m2, wav, dip);
		sf_floatwrite(dip[0], n1*n2, out);
	}

	free(wav[0][0][0]);
	free(wav[0][0]);
	free(wav[0]);
	free(wav);
	free(dip[0]);
	free(dip);
	return 0;
}



