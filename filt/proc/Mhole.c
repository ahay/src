#include <rsf.h>

int main (int argc, char* argv[])
{
    int n1, n2, n3, i1, i2, i3;
    float *pp, x,y,u,v;
    int *maskout;
    sf_file in, out, mask;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    mask = sf_output("maskout");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    sf_settype(mask,SF_INT);

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    pp = sf_floatalloc(n1);
    maskout = sf_intalloc(n1);

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) { 	
	    sf_read (pp,sizeof(float),n1,in);

	    for (i1=0; i1 < n1; i1++) { 
		x = ((float) i1)/n1 - 0.5;
		y = ((float) i2)/n2 - 0.3;
		u =  x+y;
		v = (x-y)/2.;
		if (u*u + v*v < 0.15 ) {
		    pp[i1] = 0.;
		    maskout[i1] = 0;
		} else {
		    maskout[i1] = 1;
		}
	    }
	    	    
	    sf_write (pp,sizeof(float),n1,out);
	    sf_write (maskout,sizeof(int),n1,mask);
	}
    }

    exit (0);
}


