#include <rsf.h>

#include "predict.h"
#include "banded.h"
#include "bandpass.h"
#include "prf.h"
#include "pwd.h"

static int n1, n2, nb;
static bands slv;
static float *diag, **offd, eps;
static pwd w;

void predict_init (int nx, int ny, float e, int verb)
{
    const int nw=1;

    n1 = nx;
    n2 = ny;
    nb = 2*nw;

    eps = e;

    slv = banded_init (n1, nb);
    diag = sf_floatalloc (n1);
    offd = sf_floatalloc2 (n1,nb);

    bandpass_init();
    w = pwd_init (n1, nw);
}

void predict_close (void)
{
    banded_close (slv);
    free (diag);
    free (*offd);
    free (offd);
    bandpass_close ();
    pwd_close (w);
}

void predict_flat(float** d, float** m, float** pp)
{
    int i1, i2, k2;
    float *trace, a;

    for (i2=n2-1; i2 >= 0; i2--) {
        for (i1=0; i1 < n1; i1++) {
            m[i2][i1] = d[i2][i1];
        }

        if (i2 == 0) {
            bandpass (n1, m[0]);
            return;
        }

        a = prf_burg(n1, m[i2-1]);
        prf_define (n1, a, eps, diag, offd);
        pwd_define (false, w, pp[i2-1], diag, offd);

        banded_define (slv, diag, offd);

        for (k2=n2-1; k2 >= i2; k2--) {
            trace = m[k2];

            pwd_set (w, trace, trace, diag);

            banded_solve (slv, trace);
            bandpass (n1, trace);
        }
    }
}
