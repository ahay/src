#include <valarray>

#include <rsf.hh>

int main(int argc, char* argv[])
{
    sf_init(argc,argv);

    iRSF par(0), in;
    oRSF out;

    int n1, n2;

    in.get("n1",n1);
    n2=in.size(1);

    float clip;

    par.get("clip",clip);

    std::valarray<float> trace(n1);

    for (int i2=0; i2 < n2; i2++) {
	in >> trace;

	for (int i1=0; i1 < n1; i1++) {
	    if      (trace[i1] >  clip) trace[i1]=clip;
	    else if (trace[i1] < -clip) trace[i1]=-clip;
	}

	out << trace;
    }

    exit(0);
}
