#include "rkstream.hh"
#include "gridpp_top.hh"
#include "mpigridpp.hh"
extern "C" {
#include "iwave_fopen.h"
}
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

    string fname="testsrc/mpitest1/dtestgrid.rsf";

    ofstream str;
    RVL::makeRankStream(str,rk,"testsrc/mpitest1");
    
    if (rk==0) {

      cout<<"MPIGRIDPP Unit Test 1"<<endl;
      cout<<"create MPIGridSpace, construct Vector in space"<<endl;
      cout<<"assign to perm old file=testdata.rsf"<<endl;
      cout<<"all data samples in file = 17.0, compute max and min"<<endl;
    }

    create_hfile(fname,1);
#ifdef IWAVE_USE_MPI

    MPIGridSpace<double> sp(fname,MPI_COMM_WORLD,str);
#else
    GridSpace<double> sp(fname,"notype",str);
#endif 
    Vector<double> v(sp);

    AssignFilename af(fname);
    v.eval(af);

    RVLAssignConst<double> ac(17.0);
    v.eval(ac);
    
    RVLMax<double> mx;
    MPISerialFunctionObjectRedn<double,double> mpimx(mx);
    v.eval(mpimx);
      
    RVLMin<double> mn;
    MPISerialFunctionObjectRedn<double,double> mpimn(mn);
    v.eval(mpimn);
    
    if (rk==0) cout<<"mpitest1 result: max="<<mpimx.getValue()<<" min="<<mpimn.getValue()<<endl;
    
    str.close();

#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif

    if (rk==0) iwave_fdestroy();

    return(0);

  }
  catch (RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
    exit(1);
  }
}
