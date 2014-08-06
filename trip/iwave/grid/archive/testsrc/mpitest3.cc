#include "rkstream.hh"
#include "gridpp_top.hh"
#include "mpigridpp.hh"
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

    ofstream str;
    makeRankStream(str,rk,"testsrc/mpitest3");

    /* write header file */
    string fname="testsrc/mpitest3/testgrid.rsf";

    if (rk==0) {

      cout<<"MPI GRIDPP Unit Test 3"<<endl;
      cout<<"Create GridSpace, two vectors in it, assign 1.0 to all samples,"<<endl;
      cout<<"print out inner product - should be "<<dot()<<endl;
      create_hfile(fname,0);
      
    }

#ifdef IWAVE_USE_MPI
    MPIGridSpace<float> sp(fname,MPI_COMM_WORLD,str);
#else
    GridSpace<float> sp(fname,"notype",str);
#endif    
    Vector<float> x(sp);
    Vector<float> y(sp);
      
    RVLAssignConst<float>  ac(1.0);
    x.eval(ac);
    y.eval(ac);
      
    if (rk==0) cout<<"mpitest3 result: inner product = "<<y.inner(x)<<endl;
    
    str.close();

    if (rk==0) iwave_fdestroy();

#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif

    return(0);

  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
