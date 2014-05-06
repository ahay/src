/* Test clapack (real and complex linear equations solver, cgesv_, dgesv_ */
#include<rsf.h>
#include"solver.h"

/*For the complex linear equation:*/
/*
CA=
3+3i 1+2i 2+1i
1+2i 3+3i 1+2i
2+1i 1+2i 3+3i
*/
/*
cx=
1+1i
2+2i
3+3i
*/
/*
cb=
1+21i
-4+24i
-1+27i
*/

int main(void)
{
/* 3x3 matrix A
     * 76 25 11
     * 27 89 51
     * 18 60 32
     */
/*************************************************************************************************/
/*****************          Real linear equation solving using dgesv_.           *****************/
/*************************************************************************************************/
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
/*************************************************************************************************/
/*****************    Real linear equation solving using gauss elimination.      *****************/
/*************************************************************************************************/
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
/*************************************************************************************************/
/***************** Complex linear equation solving using dgesv_ (kiss_fft_cpx).  *****************/
/*************************************************************************************************/
    int j;
    kiss_fft_cpx **ca, *cb;

    ca=(kiss_fft_cpx**) sf_complexalloc2(N,N);/*lhs matrix*/
    cb=(kiss_fft_cpx*) sf_complexalloc(N);   /*rhs and solution*/


    ca[0][0].r=3; ca[0][0].i=3;    ca[0][1].r=1; ca[0][1].i=2;    ca[0][2].r=2; ca[0][2].i=1;
    ca[1][0].r=1; ca[1][0].i=2;    ca[1][1].r=3; ca[1][1].i=3;    ca[1][2].r=1; ca[1][2].i=2;
    ca[2][0].r=2; ca[2][0].i=1;    ca[2][1].r=1; ca[2][1].i=2;    ca[2][2].r=3; ca[2][2].i=3;
    cb[0].r=1; cb[0].i=21;    cb[1].r=-4; cb[1].i=24;    cb[2].r=-1; cb[2].i=27;  
    for(j=0;j<3;j++)
      {  printf("ca [%d][0]=%.3f+i%.3f, ca [%d][1]=%.3f+i%.3f, ca [%d][2]=%.3f+i%.3f",j,ca[j][0].r,ca[j][0].i,j,ca[j][1].r,ca[j][1].i,j,ca[j][2].r,ca[j][2].i); printf("\n");}

	printf("cb [0]=%.3f+i%.3f, cb [1]=%.3f+i%.3f, cb [2]=%.3f+i%.3f",cb[0].r,cb[0].i,cb[1].r,cb[1].i,cb[2].r,cb[2].i);printf("\n");
#ifdef HAVE_MKL
    cgesv_(&N, &nrhs, (MKL_Complex8 *) ca[0], &lda, ipiv, (MKL_Complex8 *) cb, &ldb, &info);
#else
    cgesv_(&N, &nrhs, (sf_complex*) ca[0], &lda, ipiv, (sf_complex*) cb, &ldb, &info);
#endif

    printf("cx [0]=%.3f+i%.3f, cx [1]=%.3f+i%.3f, cx [2]=%.3f+i%.3f",cb[0].r,cb[0].i,cb[1].r,cb[1].i,cb[2].r,cb[2].i);printf("\n");
/*************************************************************************************************/
/***************** Complex linear equation solving using dgesv_ (kiss_fft_cpx).  *****************/
/*************************************************************************************************/

    sf_complex **ca1, *cb1;

    ca1=sf_complexalloc2(N,N);/*lhs matrix*/
    cb1=sf_complexalloc(N);   /*rhs and solution*/

    ca1[0][0]=3+3*I;  ca1[0][1]=1+2*I; ca1[0][2]=2+1*I;
    ca1[1][0]=1+2*I;  ca1[1][1]=3+3*I; ca1[1][2]=1+2*I;
    ca1[2][0]=2+1*I;  ca1[2][1]=1+2*I; ca1[2][2]=3+3*I;

    cb1[0]=1+21*I; cb1[1]=-4+24*I; cb1[2]=-1+27*I;

    /*for(j=0;j<3;j++)
        printf("ca[%d][0]=%.3f+i%.3f, ca[%d][1]=%.3f+i%.3f, ca[%d][2]=%.3f+i%.3f",j,crealf(ca1[j][0]),cimagf(ca1[j][0]),j,crealf(ca1[j][1]),cimagf(ca1[j][1]),j,crealf(ca1[j][2]),cimagf(ca1[j][2])); printf("\n");

	printf("cb[0]=%.3f+i%.3f, cb[1]=%.3f+i%.3f, cb[2]=%.3f+i%.3f",crealf(cb1[0]),cimagf(cb1[0]),crealf(cb1[1]),cimagf(cb1[1]),crealf(cb1[2]),cimagf(cb1[2])); printf("\n"); */
#ifdef HAVE_MKL
    cgesv_(&N, &nrhs, (MKL_Complex8 *) ca1[0], &lda, ipiv, (MKL_Complex8 *) cb1, &ldb, &info);
#else
    cgesv_(&N, &nrhs, ca1[0], &lda, ipiv, cb1, &ldb, &info);
#endif

    printf("cx1[0]=%.3f+i%.3f, cx1[1]=%.3f+i%.3f, cx1[2]=%.3f+i%.3f",crealf(cb1[0]),cimagf(cb1[0]),crealf(cb1[1]),cimagf(cb1[1]),crealf(cb1[2]),cimagf(cb1[2]));printf("\n");



    exit(0);
}
