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

    string fname="testsrc/test12/dtestgrid.rsf";

    if (rk==0) {

      ofstream str;
      makeRankStream(str,rk,"testsrc/test12");

      cout<<"GRIDPP Unit Test 12"<<endl;
      cout<<"create GridSpace, construct Vector in space"<<endl;
      cout<<"assign to perm old file=testgrid.rsf"<<endl;
      cout<<"all data samples in file = 12.0, compute max and min"<<endl;

      create_hfile(fname,1);

      GridSpace<double> sp(fname,"notype",str);
      Vector<double> v(sp);

      AssignFilename af(fname);
      v.eval(af);

      RVLAssignConst<double> ac(12.0);
      v.eval(ac);

      RVLMax<double> mx;
      v.eval(mx);
      
      RVLMin<double> mn;
      v.eval(mn);

      cout<<"test12 result: max="<<mx.getValue()<<" min="<<mn.getValue()<<endl;

      str.close();

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
