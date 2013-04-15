/* Test clapack */
#include<rsf.h>
#include"_lapack.h"
#include"solver.h"
int main(void)
{
/* 3x3 matrix A
     * 76 25 11
     * 27 89 51
     * 18 60 32
     */
    double A[9] = {76,27,18,25,89,60,11,51,32};
    double b[3] = {10, 7, 43};

    int N = 3;
    int nrhs = 1;
    int lda = 3;
    int ipiv[3];
    int ldb = 3;
    int info;
    
    dgesv_(&N, &nrhs, A, &lda, ipiv, b, &ldb, &info);

    if(info == 0) /* succeed */
	printf("The solution                           is %lf %lf %lf\n", b[0], b[1], b[2]);                         
    else
	printf("dgesv_ fails %d\n", info);

    float x[3], **AA;

    AA=sf_floatalloc2(N,N);
    AA[0][0]=76;  AA[0][1]=25; AA[0][2]=11;
    AA[1][0]=27;  AA[1][1]=89; AA[1][2]=51;
    AA[2][0]=18;  AA[2][1]=60; AA[2][2]=32;
    float bb[3] = {10, 7, 43};

    gaussel_init(N);
    gaussel_solve((float**) AA,bb,x);
    gaussel_close();

    printf("The solution (using gauss elimination) is %f %f %f\n", x[0], x[1], x[2]);

    exit(0);
}
