#include <rsf.h>

#include "spline3.h"

int main(int argc, char* argv[])
{
    int nd, n1, i1, two;
    float x, o1, d1, *table1, **table, *trace, x0, dx;
    bool reginput;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&two)) sf_error("Need n1= in input");
    if (two != 2) {
	reginput = true;
	nd = two;
	if (!sf_histfloat(in,"d1",&dx)) sf_error("Need d1= in input");
	if (!sf_histfloat(in,"o1",&x0)) sf_error("Need o1= in input");
    } else {
	reginput = false;
	if (!sf_histint(in,"n2",&nd)) sf_error ("Need n2= in input");
    }

    if (!sf_getint("n1",&n1)) sf_error("Need n1=");
    if (!sf_getfloat("d1",&d1)) sf_error("Need d1=");
    if (!sf_getfloat("o1",&o1)) sf_error("Need o1=");

    sf_putint(out,"n1",n1);
    sf_putfloat(out,"o1",o1);
    sf_putfloat(out,"d1",d1);
    sf_putint(out,"n2",1);

    trace = sf_floatalloc(n1);

    if (reginput) {
	spine3_init1(nd,x0,dx);
	table1 = sf_floatalloc(nd);

	sf_read(table1,sizeof(float),nd,in);
	spline_coeffs1(table1);

	for (i1=0; i1 < n1; i1++) {
	    x = o1 + i1*d1;
	    trace[i1] = spline_eval1(x);
	}
    } else {
	spine3_init(nd);
	table = sf_floatalloc2(2,nd);
	
	sf_read(table[0],sizeof(float),2*nd,in);
	spline_coeffs(table);

	for (i1=0; i1 < n1; i1++) {
	    x = o1 + i1*d1;
	    trace[i1] = spline_eval(x);
	}
    }

    sf_write(trace,sizeof(float),n1,out);

    exit(0);
}
