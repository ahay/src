#include <rsf.h>

#include "divide1.h"
#include "tridiagonal.h"

static int n;
static float eps, lam, *diag, *offd;
static tris slv;

void divide1_init(int n1, float eps1, float lam1)    
{
    int i;
    
    n = n1; eps = eps1; lam = lam1;
    slv = tridiagonal_init (n);
    diag = sf_floatalloc(n);
    offd = sf_floatalloc(n);
    for (i=0; i < n; i++) {
	offd[i] = -eps;
    }
}

void divide1_close (void)
{
    tridiagonal_close (slv);
    free (diag);
    free (offd);
}

void divide2 (int n2, float** num, float** den,  float** rat)
{
    int i2;

    divide1 (num[0], den[0], NULL, NULL, rat[0]);
    for (i2=1; i2 < n2; i2++) {
	divide1 (num[i2], den[i2], rat[i2-1], NULL, rat[i2]);
    }
}

void divide3 (int n3, int n2, float*** num, float*** den,  float*** rat)
{
    int i2, i3;
    
    divide2 (n2, num[0], den[0], rat[0]);
    for (i3=1; i3 < n3; i3++) {
	divide1 (num[i3][0], den[i3][0], rat[i3-1][0], NULL, rat[i3][0]);
	for (i2=1; i2 < n2; i2++) {
	    divide1 (num[i3][i2], den[i3][i2], 
		     rat[i3-1][i2], rat[i3][i2-1], rat[i3][i2]);
	}
    }
}

void divide1 (float* num, float* den, float* rat1, float* rat2, float* rat)
{
    int i;
    
    if (rat1 != NULL) {
	if (rat2 != NULL) {
	    for (i=0; i < n; i++) {
		diag[i] = 2.*eps + den[i]*den[i] + 2.*lam;
	    }
	    for (i=0; i < n; i++) {
		rat[i] = num[i]*den[i] + lam*(rat1[i] + rat2[i]);
	    }
	} else {
	    for (i=0; i < n; i++) {
		diag[i] = 2.*eps + den[i]*den[i] + lam;
	    }
	    for (i=0; i < n; i++) {
		rat[i] = num[i]*den[i] + lam*rat1[i];
	    }
	}
    } else {
	for (i=0; i < n; i++) {
	    diag[i] = 2.*eps + den[i]*den[i];
	}
	for (i=0; i < n; i++) {
	    rat[i] = num[i]*den[i];
	}
    }
    
    /* reflecting b.c. */
    diag[0] -= eps;
    diag[n-1] -= eps;
    
    tridiagonal_define (slv, diag, offd);
    tridiagonal_solve  (slv, rat);
}
