#include <rsf.h>

#include "predict.h"
#include "banded.h"
#include "pwd.h"

static int n1, n2, nb;
static bands slv;
static float *diag, **offd, eps;
static pwd w;

void predict_init (int nx, int ny, float e)
{
    const int nw=1;

    n1 = nx;
    n2 = ny;
    nb = 2*nw;

    eps = e;

    slv = banded_init (n1, nb);
    diag = sf_floatalloc (n1);
    offd = sf_floatalloc2 (n1,nb);

    w = pwd_init (n1, nw);
}

void predict_close (void)
{
    banded_close (slv);
    free (diag);
    free (*offd);
    free (offd);
    pwd_close (w);
}

void predict_flat(int i0, float** d, float** m, float** pp)
{
    int i1, i2, k2;
    float *trace;

    /* prediction from the left */
    for (i2=0; i2 <= i0; i2++) {
        for (i1=0; i1 < n1; i1++) {
            m[i2][i1] = d[i2][i1];
        }

        if (i2 == i0) break;

	for (i1=0; i1 < n1; i1++) {
	    diag[i1] = 6.*eps;
	    offd[0][i1] = -4.*eps;
	    offd[1][i1] = eps;
	}
	diag[0] = diag[n1-1] = 1.+eps;
	diag[1] = diag[n1-2] = 1.+5.*eps;
	offd[0][0] = offd[0][n1-2] = -2.*eps;

        pwd_define (true, w, pp[i2], diag, offd);
        banded_define (slv, diag, offd);

        for (k2=0; k2 <= i2; k2++) {
            trace = m[k2];

            pwd_set (w, trace, trace, diag);
            banded_solve (slv, trace);
        }
    }
    
    /* prediction from the right */
    for (i2=n2-1; i2 > i0; i2--) {
        for (i1=0; i1 < n1; i1++) {
            m[i2][i1] = d[i2][i1];
        }

	for (i1=0; i1 < n1; i1++) {
	    diag[i1] = 6.*eps;
	    offd[0][i1] = -4.*eps;
	    offd[1][i1] = eps;
	}
	diag[0] = diag[n1-1] = 1.+eps;
	diag[1] = diag[n1-2] = 1.+5.*eps;
	offd[0][0] = offd[0][n1-2] = -2.*eps;

        pwd_define (false, w, pp[i2-1], diag, offd);
        banded_define (slv, diag, offd);

        for (k2=n2-1; k2 >= i2; k2--) {
            trace = m[k2];

            pwd_set (w, trace, trace, diag);
            banded_solve (slv, trace);
        }
    }
}
