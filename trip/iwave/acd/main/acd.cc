#include "acd_defn.hh"
#include "istate.hh"
#include "acd_selfdoc.h"

using TSOpt::IWaveApply;
using RVL::RVLException;
int xargc;
char **xargv;

IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"csq",    0, true,  true },
  {"data",   1, false, true },
  {"source", 1, true,  false},
  {"movie",  1, false, false},
  {"init",   1, true,  false},
  {"uc_in",  1, true, true},
  {"up_in",  2, true, true},
  {"uc_out",   1, false,  true},
  {"up_out",   2, false,  true},
  {"",       0, false, false}
};

int main(int argc, char ** argv) {

  try {

    if (retrieveGlobalRank()==0 && argc<2) {
      pagedoc();
      exit(0);
    }

#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);    
#endif

    IWaveApply(argc,argv);

#ifdef IWAVE_USE_MPI
    MPI_Finalize();
#endif
  }
  catch (RVLException & e) {
    e.write(cerr);
#ifdef IWAVE_USE_MPI
    MPI_Abort(MPI_COMM_WORLD,0);
#endif
    exit(1);
  }
  
}
