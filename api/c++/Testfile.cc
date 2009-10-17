#include <valarray>

#include "rsf.hh"

int main(int argc, char* argv[])
{
    sf_init(argc,argv);

    iRSF par(0);
    iRSF in;
    oRSF out;
    int n1;
    
    in.get("n1",n1);
    if (in.type() != SF_INT)
	sf_error("Need int type.");

    std::valarray<int> trace(n1);

    in >> trace;

    out.put("n2",5);

    for (int i=0; i < 5; i++) {
	out << trace;
    }

    exit (0);
}

// 	$Id: Testfile.cc 982 2005-01-30 23:38:22Z shan $	
