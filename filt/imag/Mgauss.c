#include <math.h>

#include <rsf.h>

int main(int argc, char** argv)
{
    int n1, n2, i1, i2;
    float* vint, o1, d1, o2, d2, a, r, c1, c2, x1, x2;
    sf_file in, out;

    sf_init (argc, argv);
    in = sf_input("in");
    out = sf_output("out");
    
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint (in,"n1",&n1)) sf_error ("No n1 in input");
    if (!sf_histint (in,"n2",&n2)) sf_error ("No n2 in input");

    if (!sf_histfloat (in,"o1",&o1)) sf_error ("No o1 in input");
    if (!sf_histfloat (in,"d1",&d1)) sf_error ("No d1 in input");
    if (!sf_histfloat (in,"o2",&o2)) sf_error ("No o2 in input");
    if (!sf_histfloat (in,"d2",&d2)) sf_error ("No d2 in input");

    if (!sf_getfloat ("a",&a)) a = 1.;
    if (!sf_getfloat ("r",&r)) r = 25.;
    r = 1./(r*r);
    if (!sf_getfloat ("c1",&c1)) c1 = o1+(n1/2-1)*d1;
    if (!sf_getfloat ("c2",&c2)) c2 = o2+(n2/2-1)*d2;

    vint = sf_floatalloc(n1);
    for (i2=0; i2 < n2; i2++) {
	x2 = o2+i2*d2 - c2;
	x2 *= x2;
	sf_read (vint,sizeof(float),n1,in);
	for (i1=0; i1 < n1; i1++) {
	    x1 = o1+i1*d1 - c1;
	    x1 = (x1*x1 + x2)*r;
	    vint[i1] += a * exp(-x1);
	}
	sf_write (vint,sizeof(float),n1,out);
    }
    
    exit (0);
}
