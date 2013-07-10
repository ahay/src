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

      system("sunull nt=101 ntr=11 dt=0.002|sushw key=sx a=1000|sushw key=gx a=2000 b=-100 > testsrc/test8/hdr.su");
      
      string hdr="testsrc/test8/hdr.su";
      SEGYSpace sp(hdr);
      
      Vector<float> x(sp);

      printf("TEST 8: file system after vector creation, before filename assignment");
      iwave_fprintall(stdout);
      
      string fname="testsrc/test8/test8.su";
      AssignFilename f(fname);

      x.eval(f);

      printf("TEST 8: file system after filename assignment, before FO eval");
      iwave_fprintall(stdout);

      RVLAssignConst<float>  ac(8.0);
      x.eval(ac);
      
      printf("TEST 8: file system after FO eval");
      iwave_fprintall(stdout);

      RVLMax<float> mx;
      x.eval(mx);
      
      RVLMin<float> mn;
      x.eval(mn);
      
      cout<<"test8 result: max="<<mx.getValue()<<" min="<<mn.getValue()<<endl;
      cout<<"(should be 8, 8, and testsrc/test8/test8.su should be present)\n";

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
