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

      string cmd = "sunull nt=101 ntr=11 dt=0.002|sushw key=sx a=1000|sushw key=gx a=2000 b=-100 > testsrc/test3/hdr.su";
      system(cmd.c_str());
      
      string hdr="testsrc/test3/hdr.su";
      
      SEGYSpace sp(hdr);
      
      sp.write(cout);
      
      iwave_fdestroy();
    }

#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
    return(0);
  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
