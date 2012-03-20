/* Frequency domain Radon transform. */

#include <rsf.h>

#include "fart.h"


int main(int argc, char* argv[])
{
	int n[SF_MAX_DIM], dim;
	int n1, n2, n3;
	int i2, i3;
	float d1, d2, o1, o2;
	sf_file in, out;
	bool inv;
	int np, curv;
	float op, dp, ref, maxoff;

	sf_complex **ibuf, **obuf;

	sf_init(argc, argv);

	in  = sf_input("in");	/* f-x domain input */
	out = sf_output("out");	/* f-p domain output */

	if (SF_COMPLEX != sf_gettype(in)) 
		sf_error("FX input need complex type");

	dim = sf_filedims(in,n);
	if(dim<2) sf_error("input dim should >= 2");
	n1 = n[0];
	n2 = n[1];
	n3 = n[2];
	for(i2=3;i2<dim;i2++) n3 *= n[i2];

	if (!sf_histfloat(in,"o1",&o1)) sf_error("No o2= in input");
	if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
	if (!sf_histfloat(in,"o2",&o2)) sf_error("No o2= in input");
	if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");

	if (!sf_getint("curv",&curv)) curv=0;
	/* 0: linear; 1:parabolic */
	if (!sf_getbool("inv",&inv)) inv=false;
	/* if y, perform inverse operation */

	if (!sf_getint("no2",&np)) np=0; 
	/* stepout number */
	if (!sf_getfloat("oo2",&op)) op=0; 
	/* first stepout (moveout at "ref") */
	if (!sf_getfloat("do2",&dp)) dp=0; 
	/* stepout interval */

	maxoff = fabs((inv?(op+(np-1)*dp):(o2+d2*(n2-1))));
	if (!sf_getfloat("ref",&ref)) ref = maxoff; 
	/* reference offset */


	if(np <= 0) sf_error("no2 should larger than 0");
	sf_putint  (out, "n2", np);
	sf_putfloat(out, "o2", op);
	sf_putfloat(out, "d2", dp);

	if(inv)
	{
		ibuf = sf_complexalloc2(n1, np);
		obuf = sf_complexalloc2(n1, n2);
		sf_fart_init(curv,
			o2, d2, n2,				// stepout
			op/ref, dp/ref, np, 	// offset
			o1, d1, n1);			// freq
	}else{
		ibuf = sf_complexalloc2(n1, n2);
		obuf = sf_complexalloc2(n1, np);
		sf_fart_init(curv,
			op, dp, np, 
			o2/ref, d2/ref, n2,
			o1, d1, n1);
	}

	
	for(i3=0;i3<n3;i3++)
	{
		sf_complexread(ibuf[0], n1*n2, in);
		if(inv)	sf_ifart(ibuf, obuf);
		else	sf_fart(ibuf, obuf);
		sf_complexwrite(obuf[0], n1*np, out);
	}

	sf_fart_release();

	free(ibuf[0]);
	free(ibuf);
	free(obuf[0]);
	free(obuf);
	return 0;
}

