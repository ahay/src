#include "mpiutils.hh"
#ifdef IWAVE_USE_MPI

namespace RVL {

  /** convenient helper function to sanity-test rank settings */

  bool MPIRVL_SetRank(int & idx, int ich, MPI_Comm comm) {
    int np;
    MPI_Comm_size(comm,&np);
    if (ich < 0 || ich > np-1) {
      cerr<<"Error: MPIRVL_SetRank helper function"<<endl;
      cerr<<"assigned destination "<<ich
	  <<" out of process range [0,"<<np-1<<"]"<<cerr;
      return false;
    }
    idx=ich;;
    return true;
  }

}

#else

namespace RVL {

// dummy to avoid archiver complaints about no symbols

  void MPIRVL_foo() {}

}

#endif

