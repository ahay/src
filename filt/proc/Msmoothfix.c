#include <rsf.h>

int main(int argc, char* argv[])
{ 
    int i1, n1, i2, n2, i, nf, *fix1, *fix2;
    float **dat, **old, *fix;
    sf_file in, out;

    sf_init (argc, argv);
    in = sf_input("in");
    out = sf_output("out");

    if(!sf_histint(in,"n1",&n1)) sf_error ("No n1= in input");
    n2 = sf_leftsize(in,1);
    
    if (!sf_getint("nfix",&nf)) nf=0;

    if (nf > 0) {
	fix1 = sf_intalloc(nf);
	fix2 = sf_intalloc(nf);
	fix = sf_floatalloc(nf);

	if (!sf_getints("fix1",fix1,nf)) sf_error("Need fix1=");
	if (!sf_getints("fix2",fix2,nf)) sf_error("Need fix2=");
	if (!sf_getfloats("fix",fix,nf)) sf_error("Need fix=");
    }

    dat = sf_floatalloc2 (n1,n2);
    old = sf_floatalloc2 (n1,n2);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    dat[i2][i1] = 0.;
	}
    }
 
    for (i=0; i < nf; i++) {
	i2 = fix2[i];
	i1 = fix1[i];
	if (i1 < 0 || i1 > n1-1 || i2 < 0 || i2 > n2-1) 
	    sf_error("(%d,%d) is out of range (%d,%d)",i1,i2,n1,n2);

	dat[i2][i1] = fix[i];
    }

    sf_read (old[0],sizeof(float),n1*n2,in);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    if (dat[i2][i1] != 0.) dat[i2][i1] -= old[i2][i1];
	}
    }

    sf_write(dat[0],sizeof(float),n1*n2,out);

    exit (0);
}
