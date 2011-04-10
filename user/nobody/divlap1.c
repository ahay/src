#include <rsf.h>

#include "divlap1.h"
#include "banded.h"

/* concrete data type */
struct Divlap {
    int n;
    float eps, lam;
    bands slv;
    float *diag, *offd[2];
};

divlap divlap1_init(int n, float eps, float lam)    
{
    divlap div;
    int i;

    div = (divlap) sf_alloc(1,sizeof(*div));
    div->n = n; div->eps = eps; div->lam = lam;
    div->slv = banded_init (n,2);
    div->diag = sf_floatalloc(n);
    div->offd[0] = sf_floatalloc(n);
    div->offd[1] = sf_floatalloc(n);
    for (i=0; i < n; i++) {
	div->offd[0][i] = -4.*eps;
	div->offd[1][i] = eps;
    }

    return div;
}

void divlap1_close (divlap div)
{
    banded_close (div->slv);
    free (div->diag);
    free (div->offd[0]);
    free (div->offd[1]);
    free (div);
}

void divlap2 (divlap div, int n2, 
	      float** num, float** den, float** ref, float** rat)
{
    int i2;
    
    if (ref != NULL) {
	divlap1 (div, num[0], den[0], ref[0], NULL, NULL, rat[0]);
	for (i2=1; i2 < n2; i2++) {
	    divlap1 (div, num[i2], den[i2], ref[i2], 
		     rat[i2-1], NULL, rat[i2]);
	}
    } else {
	divlap1 (div, num[0], den[0], NULL, NULL, NULL, rat[0]);
	for (i2=1; i2 < n2; i2++) {
	    divlap1 (div, num[i2], den[i2], NULL, 
		     rat[i2-1], NULL, rat[i2]);
	}
    }
}

void divlap3 (divlap div, int n3, int n2, 
	      float*** num, float*** den, float*** ref, float*** rat)
{
    int i2, i3;
  
    divlap2 (div, n2, num[0], den[0], ref[0], rat[0]);
    for (i3=1; i3 < n3; i3++) {
	divlap1 (div, num[i3][0], den[i3][0], ref[i3][0], 
		 rat[i3-1][0], NULL, rat[i3][0]);
	for (i2=1; i2 < n2; i2++) {
	    divlap1 (div, num[i3][i2], den[i3][i2], ref[i3][i2], 
		     rat[i3-1][i2], rat[i3][i2-1], rat[i3][i2]);
	}
    }
}

void divlap1 (divlap div, float* num, float* den, float* ref, 
	      float* rat1, float* rat2, float* rat)
{
    int i;
/*    const float a=1.; */

    if (rat1 != NULL) {
	if (rat2 != NULL) {
	    for (i=0; i < div->n; i++) {
		div->diag[i] = 6.*(div->eps) + den[i] + 2.*(div->lam);
	    }
	    for (i=0; i < div->n; i++) {
		rat[i] = num[i] + (div->lam)*(rat1[i] + rat2[i]);
	    }
	} else {
	    for (i=0; i < div->n; i++) {
		div->diag[i] = 6.*(div->eps) + den[i] + (div->lam);
	    }
	    for (i=0; i < div->n; i++) {
		rat[i] = num[i] + (div->lam)*rat1[i];
	    }
	}
    } else {
	for (i=0; i < div->n; i++) {
	    div->diag[i] = 6.*(div->eps) + den[i];
	}
	for (i=0; i < div->n; i++) {
	    rat[i] = num[i];
	}
    }
    
    /* absorbing b.c.
    div->diag[0] += (a-6.)*(div->eps);
    div->diag[1] += (a-2.)*(div->eps);
    div->diag[div->n-2] += (a-2.)*(div->eps);
    div->diag[div->n-1] += (a-6.)*(div->eps);
    div->offd[0][0] = -(1.+a)*(div->eps);
    div->offd[0][div->n-2] = -(1.+a)*(div->eps);
    */

    if (ref != NULL) {
	rat[0] += (div->eps)*(ref[0]-2.*ref[1]+ref[2]);
	rat[1] += (div->eps)*(5.*ref[1]-2.*ref[0]-4.*ref[2]+ref[3]);
	for (i=2; i < div->n-2; i++) {
	    rat[i] += (div->eps)*
		(6.*ref[i]-4.*ref[i-1]-4.*ref[i+1]+ref[i-2]+ref[i+2]);
	}
	rat[div->n-2] += (div->eps)*
	    (5.*ref[div->n-2]-2.*ref[div->n-1]-4.*ref[div->n-3]+
	     ref[div->n-4]);
	rat[div->n-1] += (div->eps)*
	    (ref[div->n-1]-2.*ref[div->n-2]+ref[div->n-3]);
    }

    banded_define (div->slv, div->diag, div->offd);
    banded_solve  (div->slv, rat);
}

/* 	$Id$	 */

