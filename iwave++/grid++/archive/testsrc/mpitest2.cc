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

    ofstream str;
    makeRankStream(str,rk,"testsrc/mpitest2");

    string fname="testsrc/mpitest2/dtestgrid.rsf";

    if (rk==0) {

      cout<<"MPI GRIDPP Unit Test 2"<<endl;
      cout<<"create GridSpace, construct 2 Vectors in space"<<endl;
      cout<<"assign const 1, 2 to these, execute axpy with a=1"<<endl;
      cout<<"first vector should contain all samples = 3.0"<<endl;
      cout<<"compute max and min"<<endl;

      create_hfile(fname,1);
 
    }
    
#ifdef IWAVE_USE_MPI
    MPIGridSpace<double> sp(fname,MPI_COMM_WORLD,str);
#else
    GridSpace<double> sp(fname,"notype",str);
#endif

    Vector<double> v1(sp);
    Vector<double> v2(sp);

    RVLAssignConst<double> ac1(1.0);
    RVLAssignConst<double> ac2(2.0);
    
    v1.eval(ac1);
    v2.eval(ac2);
    
    v1.linComb(1.0,v2);
    
    RVLMax<double> mx;
    MPISerialFunctionObjectRedn<double,double> mpimx(mx);
    v1.eval(mpimx);
      
    RVLMin<double> mn;
    MPISerialFunctionObjectRedn<double,double> mpimn(mn);
    v1.eval(mpimn);
    
    if (rk==0) cout<<"mpitest2 result: max="<<mpimx.getValue()<<" min="<<mpimn.getValue()<<endl;
    
    str.close();

    if (rk==0) iwave_fdestroy();

#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif

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
