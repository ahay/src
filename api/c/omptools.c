#include "_defs.h"
#include "getpar.h"
#include "omptools.h"
#include "error.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/*------------------------------------------------------------*/
int omp_init()
/*< init OMP parameters >*/
{
    int ompnth=1;
    
#ifdef _OPENMP
    int ompath=1;

    /* OMP allowed threads */
    if(! sf_getint("ompnth",&ompnth)) ompnth=0;

#pragma omp parallel
    {	
	if(omp_get_thread_num()==0) {
	    ompath=omp_get_num_threads();
	    if(ompnth!=0)
		ompnth=SF_MIN(ompnth,ompath);
	    else
		ompnth=ompath;
	    sf_warning("using %d of %d threads",ompnth,ompath);
	}
    }
#endif 

    return ompnth;
}
