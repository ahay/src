// tests (1) segygen constructor, (2) segytrace dynamic copy
// constructor method of segygen - success if output trace
// has attributes of header file input to segygen constructor

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

      system("sunull nt=101 ntr=1 dt=0.002|sushw key=sx a=1000|sushw key=gx a=1100 > testsrc/test2/hdr.su");
      
      string hdr="testsrc/test2/hdr.su";
      SEGYDCF k(hdr);
      
      k.write(cout);
      
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
