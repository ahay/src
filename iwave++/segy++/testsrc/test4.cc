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
      
      system("sunull nt=101 ntr=11 dt=0.002|sushw key=sx a=1000|sushw key=gx a=2000 b=-100 > testsrc/test4/hdr.su");
      
      string hdr="testsrc/test4/hdr.su";
      SEGYDC d(hdr);
      
#ifdef FRUITCAKE
      cerr<<"on construction:\n";
      d.write(cerr);
#endif
      RVLAssignConst<float>  ac(1.0);
      vector<DataContainer const *> x(0);
      d.eval(ac,x);
#ifdef FRUITCAKE
      cerr<<"after evaluation:\n";
      d.write(cerr);
#endif
      
      RVLMax<float> mx;
      d.eval(mx,x);
      
      RVLMin<float> mn;
      d.eval(mn,x);
      
      cout<<"test4 result: max="<<mx.getValue()<<" min="<<mn.getValue()<<endl;
      cout<<"should be = 1, 1\n";
      cout<<"temp file not destroyed:\n";
      
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
