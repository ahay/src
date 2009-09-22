#include <rsf.h>

#include "omputil.h"

#ifdef _OPENMP
#include <omp.h>
#endif


/*------------------------------------------------------------*/
int omp_init()
/*< init OMP parameters >*/
{
    int ompnth;
    int ompchunk;
    
#ifdef _OPENMP
    int ompath;
#endif

    /* OMP data chunk size */
    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;
    
#ifdef _OPENMP
    /* OMP available threads */
    if(! sf_getint("ompnth",  &ompnth))     ompnth=0;
#pragma omp parallel
    ompath=omp_get_num_threads();
    if(ompnth<1) ompnth=ompath;
    omp_set_num_threads(ompnth);
    sf_warning("using %d threads of a total of %d",ompnth,ompath);
#else
    ompnth=0;
#endif
    
    return ompnth;
}
