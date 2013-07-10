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

    if (rk==0) {

      ofstream str;
      makeRankStream(str,rk,"testsrc/test11");

      cout<<"GRIDPP Unit Test 11"<<endl;
      cout<<"create GridDCF, use public build method from base class to"<<endl;
      cout<<"create DC, assign to perm old file=testdata.rsf"<<endl;
      cout<<"all data samples in file = 11.0, compute max and min"<<endl;

      string fname="testsrc/test11/dtestgrid.rsf";

      create_hfile(fname,1);

      GridDCF<double> f(fname,str);
      DataContainer * d = f.build();

      vector<DataContainer const *> x(0);
      
      AssignFilename af(fname);
      d->eval(af,x);

      RVLAssignConst<double> ac(11.0);
      d->eval(ac,x);

      RVLMax<double> mx;
      d->eval(mx,x);
      
      RVLMin<double> mn;
      d->eval(mn,x);

      cout<<"test11 result: max="<<mx.getValue()<<" min="<<mn.getValue()<<endl;

      delete d;

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
