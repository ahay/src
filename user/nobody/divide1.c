#include <rsf.h>

#include "divide1.h"
#include "tridiagonal.h"

/* concrete data type */
struct Div1 {
    int n;
    float eps, lam, *diag, *offd;
    tris slv;
};

div1 divide1_init(int n, float eps, float lam)    
{
    div1 div;
    int i;

    div = (div1) sf_alloc(1,sizeof(*div));
    
    div->n = n; div->eps = eps; div->lam = lam;
    div->slv = tridiagonal_init (n);
    div->diag = sf_floatalloc(n);
    div->offd = sf_floatalloc(n);
    for (i=0; i < n; i++) {
	div->offd[i] = -eps;
    }

    return div;
}

void divide1_close (div1 div)
{
    tridiagonal_close (div->slv);
    free (div->diag);
    free (div->offd);
    free (div);
}

void divide2 (div1 div, int n2, float** num, float** den,  float** rat)
{
    int i2;

    divide1 (div, num[0], den[0], NULL, NULL, rat[0]);
    for (i2=1; i2 < n2; i2++) {
	divide1 (div, num[i2], den[i2], rat[i2-1], NULL, rat[i2]);
    }
}

void divide3 (div1 div, int n3, int n2, 
	      float*** num, float*** den,  float*** rat)
{
    int i2, i3;
    
    divide2 (div, n2, num[0], den[0], rat[0]);
    for (i3=1; i3 < n3; i3++) {
	divide1 (div, num[i3][0], den[i3][0], rat[i3-1][0], NULL, rat[i3][0]);
	for (i2=1; i2 < n2; i2++) {
	    divide1 (div, num[i3][i2], den[i3][i2], 
		     rat[i3-1][i2], rat[i3][i2-1], rat[i3][i2]);
	}
    }
}

void divide1 (div1 div, float* num, float* den, 
	      float* rat1, float* rat2, float* rat)
{
    int i;
    
    if (rat1 != NULL) {
	if (rat2 != NULL) {
	    for (i=0; i < div->n; i++) {
		div->diag[i] = 2.*div->eps + den[i]*den[i] + 2.*div->lam;
	    }
	    for (i=0; i < div->n; i++) {
		rat[i] = num[i]*den[i] + div->lam*(rat1[i] + rat2[i]);
	    }
	} else {
	    for (i=0; i < div->n; i++) {
		div->diag[i] = 2.*div->eps + den[i]*den[i] + div->lam;
	    }
	    for (i=0; i < div->n; i++) {
		rat[i] = num[i]*den[i] + div->lam*rat1[i];
	    }
	}
    } else {
	for (i=0; i < div->n; i++) {
	    div->diag[i] = 2.*div->eps + den[i]*den[i];
	}
	for (i=0; i < div->n; i++) {
	    rat[i] = num[i]*den[i];
	}
    }
    
    /* reflecting b.c. */
    div->diag[0] -= div->eps;
    div->diag[div->n-1] -= div->eps;
    
    tridiagonal_define (div->slv, div->diag, div->offd);
    tridiagonal_solve  (div->slv, rat);
}

/* 	$Id$	 */
