#include "usempi.h"
#include "rkstream.hh"
#include "gridpp_top.hh"
extern "C" {
#include "iwave_fopen.h"
}
#include "wcreate_hfile.hh"
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

    string hname="testsrc/test10/testgrid.rsf";
    string fname="testsrc/test10/data10.rsf";

    if (rk==0) {

      ofstream str;
      makeRankStream(str,rk,"testsrc/test10");

      cout<<"GRIDPP Unit Test 10"<<endl;
      cout<<"create GridDC<double>, perm old file=data10.rsf"<<endl;
      cout<<"but use incompatible header file. Should throw exception"<<endl;

      create_hfile(hname,1);
      wcreate_hfile(fname,1);
      
      vector<DataContainer const *> x(0);
      GridDC<double> d(hname,str);

      AssignFilename af(fname);
      d.eval(af,x);
            
      RVLMax<double> mx;
      d.eval(mx,x);
      
      RVLMin<double> mn;
      d.eval(mn,x);

      cout<<"test10 result: max="<<mx.getValue()<<" min="<<mn.getValue()<<endl;

      str.close();

      iwave_fdestroy();

    }
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif

    return(0);

  }
  catch (RVLException & e) {
    e.write(cout);
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
    // note that this return is normal!
    exit(0);
  }
}
