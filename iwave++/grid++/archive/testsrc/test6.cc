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

    if (rk==0) {

      ofstream str;
      makeRankStream(str,rk,"testsrc/test6");

      cout<<"GRIDPP Unit Test 6"<<endl;
      cout<<"create GridDCF, GridDC, tmp file, unlink on delete, assign"<<endl;
      cout<<"6.0 to all data samples, compute max and min"<<endl;
      
      string hfile = "testsrc/test6/testgrid.rsf";

      create_hfile(hfile,0);
      
      vector<DataContainer const *> x(0);

      GridDC<float> d(hfile,str);

      RVLAssignConst<float>  ac(6.0);
      d.eval(ac,x);
      
      RVLMax<float> mx;
      d.eval(mx,x);
      
      RVLMin<float> mn;
      d.eval(mn,x);
      
      if (rk==0) {
	cout<<"test6 result: max="<<mx.getValue()<<" min="<<mn.getValue()<<endl;
	cout<<endl;
      }

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
