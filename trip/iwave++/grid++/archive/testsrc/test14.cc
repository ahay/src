#include "usempi.h"
#include "rkstream.hh"
#include "gridpp_top.hh"
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

    string fname="testsrc/test14/testgrid.rsf";

    if (rk==0) {

      ofstream str;
      makeRankStream(str,rk,"testsrc/test14");

      cout<<"GRIDPP Unit Test 14"<<endl;
      cout<<"create GridSpace, construct Vector in space"<<endl;
      cout<<"assign const value via FO eval"<<endl;
      cout<<"attempt to assign to perm old file=testdata.rsf"<<endl;
      cout<<"should throw exception, as file is already open"<<endl;

      create_hfile(fname,1);
 
      GridSpace<double> sp(fname,"notype",str);
      Vector<double> v(sp);

      RVLAssignConst<double> ac(1.0);
      v.eval(ac);

      AssignFilename af(fname);
      v.eval(af);

      str.close();

      iwave_fdestroy();
    }
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif

    return(0);

  }
  catch (RVLException & e) {
    // exception will contain run-dep filename, so just flag it
    //    e.write(cout);
    if (rk==0) 
      cout<<" *** exception thrown - exit\n";
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
    iwave_fdestroy();
    exit(0);
  }
}
