/**************** System includes **************/
#include <string.h>
#include <math.h>

/**************** RSF includes *****************/
#include <rsf.h> 

/**************** Local includes ***************/
#include "int1.h"
#include "interp_spline.h"
#include "prefilter.h"

int main(int argc, char* argv[])
{ 
    int i1, n1, i2, m2, n2, n, order, ng, ig;
    float *coord, **inp, *out, **oth, o1, d1, o2, d2, g0, dg, g;
    bool verb, noamp;
    sf_file in, warped, other;

    sf_init (argc, argv);
    in = sf_input("in");
    warped = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if(!sf_histint(in,"n1",&n1)) sf_error ("No n1= in input");
    if(!sf_histfloat(in,"d1",&d1)) sf_error ("No d1= in input");
    if(!sf_histfloat(in,"o1",&o1)) o1 = 0.;

    if(!sf_histint(in,"n2",&m2)) m2 = 1;

    if (!sf_getint("ng",&ng)) ng=1;
    if (!sf_getfloat("g0",&g0)) sf_error("Need g0=");
    if (!sf_getfloat("dg",&dg)) dg=g0;

    other = sf_input("other");

    if(!sf_histint(other,"n1",&n2)) sf_error ("No n1= in other");
    if(!sf_histfloat(other,"d1",&d2)) sf_error ("No d1= in other");
    if(!sf_histfloat(other,"o1",&o2)) o2 = 0.;

    sf_putint  (warped,"n1",n2);
    sf_putfloat(warped,"d1",d2);
    sf_putfloat(warped,"o1",o2);

    sf_putint  (warped,"n3",ng);
    sf_putfloat(warped,"d3",dg);
    sf_putfloat(warped,"o3",g0);

    n = n2*m2;

    if(!sf_getbool("verb",&verb)) verb = false;
    if(!sf_getbool("noamp",&noamp)) noamp = false;

    if(!sf_getint("accuracy",&order)) {
	order = 2;
    } else if (order < 1 || order > 4) {
	sf_error ("accuracy must be between 1 and 4");
    }
    order *= 2;

    coord = sf_floatalloc (n2); 
    inp =   sf_floatalloc2 (n1,m2);
    out =   sf_floatalloc (n2);
    oth =   sf_floatalloc2 (n2,m2);

    prefilter_init (order, n1, order*10);     
    for (i2=0; i2 < m2; i2++) {
	sf_read(inp[i2],sizeof(float),n1,in);
	prefilter_apply (n1, inp[i2]);
    }
    prefilter_close();

    sf_read(oth[0],sizeof(float),m2*n2,other);
    sf_fileclose(other);

    for (ig=0; ig < ng; ig++) {
	g = g0 + ig*dg;

	for (i1=0; i1 < n2; i1++) {
	    coord[i1] = (o2+i1*d2)*g;
	}

	int1_init (coord, o1, d1, n1, spline_int, order, n2);

	for (i2=0; i2 < m2; i2++) {
	    int1_lop (false,false,n1,n2,inp[i2],out);
	    for (i1=0; i1 < n2; i1++) {
		out[i1] -= oth[i2][i1];
	    }
	    sf_write(out,sizeof(float),n2,warped);
	}
    }

    exit (0);
}
