#include <rsf.h>

#include "helix.h"
#include "nhelix.h"

nfilter nallocate(int np, int nd, int *nh, int *pch) 
{
    nfilter aa;
    int ip, id;
    
    aa = (nfilter) sf_alloc(1,sizeof(*aa));

    aa->np = np;
    aa->hlx = (filter*) sf_alloc(np,sizeof(filter));

    for (ip=0; ip < np; ip++) {
	aa->hlx[ip] = allocatehelix(nh[ip]);
    }
    
    aa->pch = sf_intalloc(nd);
    for (id=0; id < nd; id++) {
	aa->pch[id] = pch[id];
    }

    aa->mis = NULL;
    return aa;
}

void ndeallocate(nfilter aa) 
{
    int ip;
    
    for (ip=0; ip < aa->np; ip++) {
	deallocatehelix(aa->hlx[ip]);
    }
    free (aa->hlx);
    free (aa->pch);
    if (NULL != aa->mis) free (aa->mis);
    free (aa);
}

