#include <rsf.h>

#include "createnhelix.h"
#include "createhelix.h"
#include "nhelix.h"
#include "nbound.h"

nfilter createnhelix(int dim, int *nd, int *center, int *gap, int *na, 
		     int *pch) 
{
    nfilter nsaa;
    filter aa;
    int n123, np, ip, *nh, i;

    aa = createhelix(dim, nd, center, gap, na);

    n123=1;
    for (i=0; i < dim; i++) {
	n123 *= nd[i];
    }
    np = pch[0];
    for (i=0; i < n123; i++) {
	if (pch[i] > np) np=pch[i];
    }

    nh = sf_intalloc(np);
    for (ip=0; ip < np; ip++) {
	nh[ip] = aa->nh;
    }
    nsaa = nallocate(np, n123, nh, pch);
    for (ip=0; ip < np; ip++) {
	for (i=0; i < aa->nh; i++) {
	    nsaa->hlx[ip]->lag[i] = aa->lag[i];
	}
	deallocatehelix(aa);
	nbound(ip, dim, nd, na, nsaa); 
    }

    return nsaa;
}

