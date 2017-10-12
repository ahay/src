// Clip the data.

#include <rsf.hh>
#include <float.h>
#include <valarray>

int main(int argc, char* argv[])
{
    // trace length, number of traces
    int n1, n2;
    float upper, lower;
    
    // Initialize RSF
    sf_init(argc,argv); 
    
    // input parameter, file
    iRSF par(0), in;
    // output file
    oRSF out;

	// check that the input is float
	if (SF_FLOAT != in.type())
	sf_error("Need float input");

	// n1 is the fastest dimension
    in.get("n1",n1);

	// leftsize gets n2*n3*n4*...
    n2=in.size(1);

    // parameter from the command line
    par.get("upper",upper,+FLT_MAX);
    par.get("lower",lower,-FLT_MAX);

	// allocate floating point array
    std::valarray<float> trace(n1);

	// loop over traces
    for (int i2=0; i2 < n2; i2++) {
	
	// read a trace (overloading operator)	
	in >> trace;

	// loop over samples
	for (int i1=0; i1 < n1; i1++) {
	    if      (trace[i1] > upper) trace[i1]=upper;
	    else if (trace[i1] < lower) trace[i1]=lower;
	}

	// write a trace
	out << trace;
    }

    exit(0);
}
