/* Simple 2-D synthetics with crossing plane waves.

Takes: > synth.rsf
*/
#include <math.h>

#include <rsf.h>

#include "triangle.h"
#include "random.h"

int main(int argc, char* argv[])
{
    int i1, i2, i3, n1, n2, n3, n, ns, p, t1, t2;
    float **pp, *s1, *s2;
    triangle tr1, tr2;
    sf_file mod;

    sf_init (argc,argv);
    mod = sf_output("out");
    sf_setformat(mod,"native_float");

    if (!sf_getint("n1",&n1)) n1=100;
    if (!sf_getint("n2",&n2)) n2=14;
    if (!sf_getint("n3",&n3)) n3=1;
    
    sf_putint(mod,"n1",n1);
    sf_putint(mod,"n2",n2);
    sf_putint(mod,"n3",n3);

    pp = sf_floatalloc2(n1,n2);

    if (!sf_getint("n",&n)) n=3;
    if (!sf_getint("p",&p)) p=3;
    ns = n1*n;
    
    s1 = sf_floatalloc(ns);
    s2 = sf_floatalloc(ns);

    random_init (1993);
    for (i1=0; i1 < ns; i1++) {
	s1[i1] = powf(random0()-0.5,(float) p);
	s2[i1] = powf(random0()-0.5,(float) p);
    }

    if (!sf_getint("t1",&t1)) t1=4;
    if (!sf_getint("t2",&t2)) t2=4;

    tr1 = triangle_init (t1, ns);
    tr2 = triangle_init (t2, ns);

    smooth2(tr1,0,1,false,s1);
    smooth2(tr2,0,1,false,s2);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    pp[i2][i1]  = s1[i1 + n/2 * n1 + i2];
	    pp[i2][i1] += s2[i1 + n/2 * n1 - i2*2];
	}
    }

    for (i3=0; i3 < n3; i3++) {
	sf_floatwrite (pp[0],n1*n2,mod);
    }
    
    exit(0);
}

