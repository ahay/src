#include <rsf.h>

#include "repeat.h"

static int n1, n2;
static sf_operator oper;

void repeat_init(int m1, int m2, sf_operator oper1)
{
    n1 = m1;
    n2 = m2;
    oper = oper1;
}

/* Causal integration */
void repeat_lop (bool adj, bool add, int nx, int ny, float *xx, float *yy)
{
    int i2;       
    
    if (nx != ny || nx != n1*n2) sf_error("%s: Wrong size",__FILE__);

    sf_adjnull (adj, add, nx, ny, xx, yy);

    for (i2=0; i2 < n2; i2++) {
	oper(adj,true,n1,n1,xx+i2*n1,yy+i2*n1);
    }
}

/* 	$Id: repeat.c,v 1.1 2004/04/05 14:38:29 fomels Exp $	 */

