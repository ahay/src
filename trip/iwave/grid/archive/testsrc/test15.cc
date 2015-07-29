#include "usempi.h"
#include "rkstream.hh"
#include "gridpp_top.hh"
extern "C" {
#include "iwave_fopen.h"
}
#include "create_hfile.hh"

using namespace RVL;
using namespace TSOpt;

extern float dot();

int main(int argc, char ** argv) {

  int rk=0;

#ifdef IWAVE_USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rk);
#endif 

  try {

    
    /* write header file */
    string fname="testsrc/test15/testgrid.rsf";

    if (rk==0) {

      ofstream str;
      makeRankStream(str,rk,"testsrc/test15");

      cout<<"GRIDPP Unit Test 15"<<endl;
      cout<<"Create GridSpace, two vectors in it, assign 1.0 to all samples,"<<endl;
      cout<<"print out inner product - should be "<<dot()<<endl;

      create_hfile(fname,0);

      GridSpace<float> sp(fname,"notype",str);
      
      Vector<float> x(sp);
      Vector<float> y(sp);
      
      RVLAssignConst<float>  ac(1.0);
      x.eval(ac);
      y.eval(ac);
      
      cout<<"test11 result: inner product = "<<y.inner(x)<<endl;
      cout<<endl;
    
      str.close();

      iwave_fdestroy();

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
