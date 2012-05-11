/* omnidirectional dip estimation  */

#include <rsf.h>
#include "odip.h"

int main(int argc, char*argv[])
{
	sf_file in, out;
	int nf, n1, n2, n3, rect[2], niter, liter;
	int i3;
	bool verb, opwd;
	float **wav, **dip;

	sf_init(argc, argv);

	in  = sf_input("in");
	out = sf_output("out");

	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

	if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
	if (!sf_histint(in, "n2", &n2)) sf_error("No n2= in input");
	n3 = sf_leftsize(in, 2);

	if (!sf_getint("nf",&nf)) nf=1;
	/* PWD filter order */
	if (!sf_getint("rect1",&rect[0])) rect[0]=1;
	/* dip smoothness on 1st axis */
	if (!sf_getint("rect2",&rect[1])) rect[1]=1;
	/* dip smoothness on 2nd axis */
	
	if (!sf_getint("niter",&niter)) niter=5;
	/* number of iterations */
	if (!sf_getint("liter",&liter)) liter=20;
	/* number of linear iterations */

	if (!sf_getbool("opwd", &opwd)) verb = true;
	/* y: opwd;  n: lpwd */
	if (!sf_getbool("verb",&verb)) verb = false;
	/* verbosity flag */

	wav = sf_floatalloc2(n1,n2);
	dip = sf_floatalloc2(n1,n2);

	/* initialize dip estimation */
	odip_init(nf, n1, n2, rect, liter, verb);
	for(i3=0; i3<n3; i3++)
	{
		sf_floatread(wav[0], n1*n2, in);
		odip(wav, dip, niter, opwd);
		sf_floatwrite(dip[0], n1*n2, out);
	}

	odip_close();
	return 0;
}

