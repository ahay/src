#include "asg_defn.hh"
#include "istate.hh"
#include "asg_selfdoc.h"

using TSOpt::IWaveApply;
using RVL::RVLException;

int xargc_;
char **xargv_;

IOKEY IWaveInfo::iwave_iokeys[]
= {
  {"bulkmod",    0, true,  true },
  {"buoyancy",   1, true,  true },
  {"source_p",   2, true,  false},
  {"data_p",     2, false, true},
  {"data_v0",    5, false, true},
  {"data_v1",    6, false, true},
  {"data_v2",    7, false, true},
  {"movie_p",    2, false, false},
  {"movie_v0",   5, false, false},
  {"movie_v1",   6, false, false},
  {"movie_v2",   7, false, false},
  {"",           0, false, false}
};

int main(int argc, char ** argv) {

  try {

#ifdef IWAVE_USE_MPI
    int ts=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&ts);    
#endif

    if (retrieveGlobalRank()==0 && argc<2) {
      pagedoc();
      exit(0);
    }

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
