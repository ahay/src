/* utilities for openmp warmup*/

#include <rsf.h>
#ifdef _OpenMP
#include <omp.h>
#endif

#ifndef _newton_h

typedef float(*funcl)(float);
/*function pointer for float -> float */

#endif

void printhreads(){
  
  #ifdef _OpenMP
  #pragma omp critical
    sf_warning("Hello from thread %d",omp_get_thread_num());
  #endif
  return;
}
