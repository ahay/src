#include <rsf.h>

#include "pick0.h"

int main (int argc, char *argv[])
{
    int i1, n1, n2, n3, i2, i3, order;
    float *dip, *pik;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    /* transpose */
    sf_putint(out,"n1",n2);
    sf_putint(out,"n2",n1);
    sf_putfloat(out,"o2",0.);
    sf_putfloat(out,"d2",1.);

    if (!sf_getint("order",&order)) order=4;
    
    dip = sf_floatalloc (n1);
    pik = sf_floatalloc (n2);

    pick0_init (n1, n2, order);

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    sf_read(dip,sizeof(float),n1,in);
	    pick0_set (i2, dip);
	}

	for (i1=0; i1 < n1; i1++) {
	    pick0_step0 ((float) i1, pik);
	    sf_write(pik,sizeof(float),n2,out);
	}
    }

    exit (0);
}
