#include <rsf.h>

#include "divide2.h"
#include "tridiagonal.h"

/* concrete data type */
struct Div2 {
    int n1, n2;
    float eps1, eps2;
    float *diag1, *diag2, *offd1, *offd2, *rat2;
    tris slv1, slv2;
};

div2 divide2_init(int n1, int n2, float eps1, float eps2)    
{
    div2 div;
    int i;

    div = (div2) sf_alloc(1,sizeof(*div));
    
    eps1 *= 2;
    eps2 *= 2;

    div->n1 = n1;
    div->n2 = n2;
    div->eps1 = eps1; 
    div->eps2 = eps2;
    div->slv1 = tridiagonal_init (n1);
    div->slv2 = tridiagonal_init (n2);
    div->diag1 = sf_floatalloc(n1);
    div->diag2 = sf_floatalloc(n2);
    div->offd1 = sf_floatalloc(n1);
    div->offd2 = sf_floatalloc(n2);
    for (i=0; i < n1; i++) {
	div->offd1[i] = -eps1;
    }
    for (i=0; i < n2; i++) {
	div->offd2[i] = -eps2;
    }
    div->rat2 = sf_floatalloc(n2);

    return div;
}

void divide2_close (div2 div)
{
    tridiagonal_close (div->slv1);
    tridiagonal_close (div->slv2);
    free (div->diag1);
    free (div->diag2);
    free (div->offd1);
    free (div->offd2);
    free (div->rat2);
    free (div);
}

/* AOS - additive splitting */
void divide2 (div2 div, float** num, float** den,  float** rat)
{
    int i1, i2, n1, n2;
    float deni;

    n1 = div->n1;
    n2 = div->n2;

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    deni = den[i2][i1];
	    div->diag1[i1] = 2.*div->eps1 + deni*deni;
	    rat[i2][i1] = 0.5*num[i2][i1]*deni;
	}
    
	/* almost reflecting b.c. */
	div->diag1[0] -= 0.999*div->eps1;
	div->diag1[n1-1] -= 0.999*div->eps1;
    
	tridiagonal_define (div->slv1, div->diag1, div->offd1);
	tridiagonal_solve  (div->slv1, rat[i2]);
    }
    for (i1=0; i1 < n1; i1++) {
	for (i2=0; i2 < n2; i2++) {
	    deni = den[i2][i1];
	    div->diag2[i2] = 2.*div->eps2 + deni*deni;
	    div->rat2[i2] = 0.5*num[i2][i1]*deni;
	}
    
	/* almost reflecting b.c. */
	div->diag2[0] -= 0.999*div->eps2;
	div->diag2[n2-1] -= 0.999*div->eps2;
    
	tridiagonal_define (div->slv2, div->diag2, div->offd2);
	tridiagonal_solve  (div->slv2, div->rat2);

	for (i2=0; i2 < n2; i2++) {
	    rat[i2][i1] += div->rat2[i2];
	}
    }
}

/* 	$Id$	 */
