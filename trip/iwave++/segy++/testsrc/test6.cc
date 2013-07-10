#include "usempi.h"
#include "functions.hh"
#include "segypp.hh"

//#define FRUITCAKE

using namespace RVL;
using namespace TSOpt;

char ** xargv;

int main(int argc, char ** argv) {

  int rk=0;

#ifdef IWAVE_USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rk);
#endif 

  try {

    if (rk==0) {

      system("sunull nt=101 ntr=11 dt=0.002|sushw key=sx a=1000|sushw key=gx a=2000 b=-100 > testsrc/test6/hdr.su");
      system("cp testsrc/test6/hdr.su testsrc/test6/test6.su");
      
      string hdr="testsrc/test6/hdr.su";
      SEGYDC d(hdr);

#ifdef FRUITCAKE
      cerr<<"on construction:\n";
      d.write(cerr);
#endif
      
      // dummy vector of DC * args
      vector<DataContainer const *> x(0);
      
      AssignFilename f("testsrc/test6/test6.su");
      d.eval(f,x);
      
#ifdef FRUITCAKE
      cerr<<"after AssignFile evaluation:\n";
      d.write(cerr);
#endif
      
      RVLAssignConst<float>  ac(6.0);
      d.eval(ac,x);
      
#ifdef FRUITCAKE
      cerr<<"after AssignConst evaluation:\n";
      d.write(cerr);
#endif
      
      RVLMax<float> mx;
      d.eval(mx,x);
      
      RVLMin<float> mn;
      d.eval(mn,x);
      
      cout<<"test6 result: max="<<mx.getValue()<<" min="<<mn.getValue()<<endl;
      cout<<"(should be 6, 6, and testsrc/hdr.su, testsrc/test6.su should exist\n";

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
