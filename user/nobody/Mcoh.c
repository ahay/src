/* 3-D coherency estimation using plane wave destruction.

Takes: < data.rsf > residual.rsf
*/

#include <rsf.h>

#include "coh.h"

int main (int argc, char *argv[])
{
    int n1,n2,n3, w1,w2,w3, m1,m2,m3, i1,i2,i3, n123, nw, nj1, nj2;
    float ***u, ***p, ***c;
    sf_file in, out, dip;
    coh ch;

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
    /* [1,2,3] accuracy */
    if (nw < 1 || nw > 3) 
	sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);
    if (!sf_getint("nj1",&nj1)) nj1=1;
    /* in-line aliasing */
    if (!sf_getint("nj2",&nj2)) nj2=1;
    /* cross-line aliasing */

    if (!sf_getint("w1",&w1)) sf_error("Need w1=");
    if (!sf_getint("w2",&w2)) sf_error("Need w2=");
    if (!sf_getint("w3",&w3)) w3=1;
    /* window size */
    if (w1 > n1) w1=n1;
    if (w2 > n2) w2=n2;
    if (w3 > n3) w3=n3;

    u = sf_floatalloc3(n1,n2,n3);
    p = sf_floatalloc3(n1,n2,n3);
    c = sf_floatalloc3(n1,n2,n3);

    /* read data */
    sf_floatread(u[0][0],n123,in);

    /* read t-x dip */
    sf_floatread(p[0][0],n123,dip);

    ch = coh_init (nw,nj1,n1,n2,n3,p);
  
    /* apply */
    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		c[i3][i2][i1] = coh1(ch,i1,i2,i3,w1,w2,w3,u);
	    }
	}
    }

    /* write t-x coherency */
    sf_floatwrite(c[0][0],n123,out);

    if (n3 > 1) { /* if 3-D input */
	/* read t-y dip */
	sf_floatread(p[0][0],n123,dip);
	
	if (nj2 != nj1) ch = coh_init (nw,nj2,n1,n2,n3,p);
  
	/* apply */
	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    c[i3][i2][i1] = coh2(ch,i1,i2,i3,w1,w2,w3,u);
		}
	    }
	}

	/* write t-y coherency */
	sf_floatwrite(c[0][0],n123,out);
    }
    

    exit (0);
}

/* 	$Id$	 */
