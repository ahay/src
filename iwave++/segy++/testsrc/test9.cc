#include "usempi.h"
#include "segypp.hh"

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

      system("sunull nt=101 ntr=11 dt=0.002|sushw key=sx a=1000|sushw key=gx a=2000 b=-100 > testsrc/test9/hdr.su");
      
      string hdr="testsrc/test9/hdr.su";
      SEGYSpace sp(hdr);

      Vector<float> x(sp);
    
      RVLAssignConst<float>  ac(9.0);
      x.eval(ac);
      
      RVLMax<float> mx;
      x.eval(mx);
      
      RVLMin<float> mn;
      x.eval(mn);
      
      cout<<"test9 result: max="<<mx.getValue()<<" min="<<mn.getValue()<<endl;
      cout<<"(should be 9, 9, temp file should be unlinked)\n";
      
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
