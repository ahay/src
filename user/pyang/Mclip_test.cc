/* Clip the data. */

#include <valarray>
#include <rsf.hh>

int main(int argc, char* argv[])
{
    sf_init(argc,argv); //! Initialize RSF
    
    iRSF par(0), in; //! input parameter, file
    oRSF out;        //! output file

    int n1, n2;      //!trace length, number of traces
    float clip;
    
    in.get("n1",n1);
    n2=in.size(1);

    par.get("clip",clip); // parameter from the command line

    std::valarray<float> trace(n1);

    for (int i2=0; i2 < n2; i2++) { // loop over traces
	in >> trace; // read a trace

	for (int i1=0; i1 < n1; i1++) { // loop over samples
	    if      (trace[i1] >  clip) trace[i1]=clip;
	    else if (trace[i1] < -clip) trace[i1]=-clip;
	}

	out << trace; // write a trace
    }

    exit(0);
}
