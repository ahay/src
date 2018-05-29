/* utilities for openmp warmup*/

#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef _ompinitialize_h

typedef float(*funcl)(float);
/*function pointer for float -> float */
/*^*/

#endif

void printthreads(){
  
  #ifdef _OPENMP
  #pragma omp critical
    sf_warning("Hello from thread %d",omp_get_thread_num());
  #endif
  return;
}
