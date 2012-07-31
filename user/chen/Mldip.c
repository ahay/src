/* dip estimation by line interpolating pwd */

#include <rsf.h>
#include "ldip.h"


int main(int argc, char*argv[])
{
	sf_file in, out, p0;
	int nf, n1, n2, n3, rect[2], niter, liter, interp;
	int i3;
	bool verb;
	float **wav, **dip, eta, **dip0;

	sf_init(argc, argv);

	in  = sf_input("in");
	out = sf_output("out");

	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

	if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
	if (!sf_histint(in, "n2", &n2)) sf_error("No n2= in input");
	n3 = sf_leftsize(in, 2);

	if (!sf_getint("rect1",&rect[0])) rect[0]=0;
	/* dip smoothness on 1st axis */
	if (!sf_getint("rect2",&rect[1])) rect[1]=0;
	/* dip smoothness on 2nd axis */
	
	if (!sf_getint("order",&nf)) nf=1;
	/* pwd order */
	if (!sf_getint("niter",&niter)) niter=5;
	/* number of iterations */
	if (!sf_getint("liter",&liter)) liter=20;
	/* number of linear iterations */
	if (!sf_getint("interp",&interp)) interp=0;
	/* interpolation method: 
	0: maxflat
	1: Lagrange 
	2: B-Spline */
	if (!sf_getfloat("eta", &eta)) eta = 1.0;
	/* steps for iteration */

	if(sf_getstring("dip0")!=NULL)
	{
		dip0 = sf_floatalloc2(n1, n2);
		p0 = sf_input("dip0");
		sf_histint(p0, "n1", &i3);
		if(i3!=n1) sf_error("p0.n1=%d and in.n1=%d",i3,n1);
		sf_histint(p0, "n2", &i3);
		if(i3!=n2) sf_error("p0.n2=%d and in.n2=%d",i3,n2);
		sf_floatread(dip0[0], n1*n2, p0);
	}else	p0=NULL;
	
	/* starting dip */
	if (!sf_getbool("verb", &verb)) verb = false;
	/* verbosity flag */

	wav = sf_floatalloc2(n1, n2);
	dip = sf_floatalloc2(n1, n2);

	

	/* initialize dip estimation */
	ldip_init(interp, nf, n1, n2, rect, liter, verb);


	for(i3=0; i3<n3; i3++)
	{
		sf_floatread(wav[0], n1*n2, in);
		if(p0!=NULL) memcpy(dip[0], dip0[0], n1*n2*sizeof(float));
		else memset(dip[0], 0, n1*n2*sizeof(float));
		ldip(wav, dip, niter, eta);
		sf_floatwrite(dip[0], n1*n2, out);
	}

	if(p0!=NULL)
	{
		free(dip0[0]);
		free(dip0);
	}
	ldip_close();
	free(dip[0]);
	free(wav[0]);
	free(dip);
	free(wav);
	return 0;
}



