#include <rsf.h>

#include "allp3.h"

int main (int argc, char *argv[])
{
    int n1,n2,n3, m1, m2, m3, n123, nw, nj1, nj2;
    float ***u1, ***u2, ***p;
    sf_file in, out, dip;
    allpass ap;

    sf_init(argc,argv);
    in = sf_input ("in");
    dip = sf_input ("dip");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in) ||
	SF_FLOAT != sf_gettype(dip)) sf_error("Need float type");

    if (!sf_histint(in,"n1",&n1)) sf_error("Need n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    if (!sf_histint(in,"n3",&n3)) n3=1;
    n123 = n1*n2*n3;

    if (!sf_histint(dip,"n1",&m1) || m1 != n1) 
	sf_error("Need n1=%d in dip",n1);
    if (1 != n2 && (!sf_histint(dip,"n2",&m2) || m2 != n2)) 
	sf_error("Need n2=%d in dip",n2);
    if (1 != n3 && (!sf_histint(dip,"n3",&m3) || m3 != n3)) 
	sf_error("Need n3=%d in dip",n3);

    /* two dips output in 3-D */
    if (n3 > 1) sf_putint(out,"n4",2); 

    if (!sf_getint("order",&nw)) nw=1;
    if (nw < 1 || nw > 3) 
	sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);
    if (!sf_getint("nj1",&nj1)) nj1=1;
    if (!sf_getint("nj2",&nj2)) nj2=1;

    u1 = sf_floatalloc3(n1,n2,n3);
    u2 = sf_floatalloc3(n1,n2,n3);
    p  = sf_floatalloc3(n1,n2,n3);

    /* read data */
    sf_read(u1[0][0],sizeof(float),n123,in);

    /* read t-x dip */
    sf_read(p[0][0],sizeof(float),n123,dip);

    ap = allpass_init (nw,nj1,n1,n2,n3,p);
  
    /* apply */
    allpass1(false, ap, u1, u2);

    /* write t-x destruction */
    sf_write(u2[0][0],sizeof(float),n123,out);

    if (n3 > 1) { /* if 3-D input */
	/* read t-y dip */
	sf_read(p[0][0],sizeof(float),n123,dip);
	
	if (nj2 != nj1) ap = allpass_init (nw,nj2,n1,n2,n3,p);
  
	/* apply */
	allpass2(false, ap, u1, u2);

	/* write t-y destruction */
	sf_write(u2[0][0],sizeof(float),n123,out);
    }
    
    exit (0);
}
