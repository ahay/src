#ifndef __RVL_MPI_UTILS
#define __RVL_MPI_UTILS

#ifdef IWAVE_USE_MPI

#include "mpi.h"
#include "std_cpp_includes.hh"

namespace RVL {


  /** convenient helper function to sanity-test rank settings */

  bool MPIRVL_SetRank(int & idx, int ich, MPI_Comm comm);

}
#endif
#endif
