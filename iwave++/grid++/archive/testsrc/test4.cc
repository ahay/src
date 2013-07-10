#include "usempi.h"
#include "rkstream.hh"
#include "gridpp_top.hh"
#include "create_hfile.hh"

using namespace RVL;
using namespace TSOpt;

int main(int argc, char ** argv) {

  int rk=0;

#ifdef IWAVE_USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rk);
#endif 

  try {

    string fname="testsrc/test4/testgrid.rsf";

    ofstream str;
    makeRankStream(str,rk,"testsrc/test4");

    if (rk==0) {

      cout<<"GRIDPP Unit Test 4"<<endl;
      cout<<"construct GridDCF, write method\n";

      create_hfile(fname,0);

      GridDCF<float> f(fname,str);
      
      f.write(cout);
      cout<<endl;

    }
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif

    return(0);

  }
  catch (RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
}
