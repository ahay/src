#include "prf.h"

float prf_burg(int n, float* trace)
{
    int i;
    float a, avto, cros;

    avto = trace[n-1]*trace[n-1];
    a = avto + trace[0]*trace[0];
    cros = 0.;
    for (i=0; i < n-1; i++) {
        avto += trace[i]*trace[i];
        cros += trace[i]*trace[i+1];
    }
    if (avto == 0.) {
        a=-2.;
    } else {
        a = - (a + 2.*cros)/avto;
    }

    return a;
}

float prf_burg2(int n1, int n2, float** trace)
{
    int i1, i2;
    float a;
    double avto, cros;

    a = avto = cros = 0.;

    for (i2=0; i2 < n2; i2++) {
        avto += trace[i2][n1-1]*trace[i2][n1-1];
        a += (trace[i2][n1-1]*trace[i2][n1-1] + trace[i2][0]*trace[i2][0]);
        for (i1=0; i1 < n1-1; i1++) {
            avto += trace[i2][i1]*trace[i2][i1];
            cros += trace[i2][i1]*trace[i2][i1+1];
        }
    }

    if (avto == 0.) {
        a=-2.;
    } else {
        a = - (a + 2.*cros)/avto;
    }

    return a;
}

void prf_define (int n, float a, float eps, float* diag, float** offd)
{
    int i;

    for (i=0; i < n-2; i++) {
        diag[i+1] = (2. + a*a)*eps;
        offd[0][i] = (2.*a)*eps;
        offd[1][i] = eps;
    }
    offd[0][0]   = (1.+2.*a)*eps;
    offd[0][n-2] = (1.+2.*a)*eps;
    diag[0]    = (1.+(1.+a)*(1.+a))*eps;
    diag[n-1]  = (1.+(1.+a)*(1.+a))*eps;
}
