// tests basic container package for segy's - CP<float,segy>=segytrace

#include "usempi.h"
#include "segypp.hh"
#include "segyfun.hh"

using namespace TSOpt;

int main(int argc, char ** argv) {

  int rk=0;

#ifdef IWAVE_USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rk);
#endif 

  try {

    if (rk==0) {
      
      segytrace a;
      
      segy m;
      
      m.ns=101;
      m.dt=4000;
      m.delrt=200;
      
      cout<<"initialzation success = "<<a.initialize(m)<<endl;
      
      cout<<"function: data size = "<<getDataSize<segy>(m)<<endl;
      
      cout<<"class:    data size = "<<a.getSize()<<endl;
   
    }
#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
    return 0;
  }
  catch (RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
}
