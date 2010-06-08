#ifndef _BLAS_H_
#define _BLAS_H_

//#include "bltypes.h"
#include <complex>
typedef std::complex<float> cpx8;
typedef std::complex<double> cpx16;

extern "C"
{
  float sasum_(int *n,float *x,int *incx);
  void  saxpy_(int *n,float *alpha,float *x,int *incx,float *y,int *incy);
  void  saxpyi_(int *nz,float *a,float *x,int *indx,float *y);
  float scasum_(int *n,cpx8 *x,int *incx); 
  float scnrm2_(int *n,cpx8 *x,int *incx); 
  void  scopy_(int *n,float *x,int *incx,float *y,int *incy);
  float sdot_(int *n,float *x,int *incx,float *y,int *incy);
  float sdoti_(int *nz,float *x,int *indx,float *y);
  void  sgthr_(int *nz,float *y,float *x,int *indx);
  void  sgthrz_(int *nz,float *y,float *x,int *indx);
  float snrm2_(int *n,float *x,int *incx);
  void  srot_(int *n,float *x,int *incx,float *y,int *incy,float *c,float *s);
  void  srotg_(float *a,float *b,float *c,float *s);
  void  sroti_(int *nz,float *x,int *indx,float *y,float *c,float *s);
  void  srotm_(int *n,float *x,int *incx,float *y,int *incy,float *param);
  void  srotmg_(float *d1,float *d2,float *x1,float *y1,float *param);
  void  sscal_(int *n,float *a,float *x,int *incx);
  void  ssctr_(int *nz,float *x,int *indx,float *y);
  void  sswap_(int *n,float *x,int *incx,float *y,int *incy);
  int   isamax_(int *n,float *x,int *incx);
  int   isamin_(int *n,float *x,int *incx);
  
  void caxpy_(int *n,cpx8 *alpha,cpx8 *x,int *incx,cpx8 *y,int *incy); 
  void caxpyi_(int *nz,cpx8 *a,cpx8 *x,int *indx,cpx8 *y); 
  void ccopy_(int *n,cpx8 *x,int *incx,cpx8 *y,int *incy); 
  void cdotc_(cpx8 *pres,int *n,cpx8 *x,int *incx,cpx8 *y,int *incy); 
  void cdotci_(cpx8 *pres,int *nz,cpx8 *x,int *indx,cpx8 *y); 
  void cdotu_(cpx8 *pres,int *n,cpx8 *x,int *incx,cpx8 *y,int *incy); 
  void cdotui_(cpx8 *pres,int *nz,cpx8 *x,int *indx,cpx8 *y); 
  void cgthr_(int *nz,cpx8 *y,cpx8 *x,int *indx); 
  void cgthrz_(int *nz,cpx8 *y,cpx8 *x,int *indx); 
  void crotg_(cpx8 *a,cpx8 *b,float *c,cpx8 *s); 
  void cscal_(int *n,cpx8 *a,cpx8 *x,int *incx); 
  void csctr_(int *nz,cpx8 *x,int *indx,cpx8 *y); 
  void csrot_(int *n,cpx8 *x,int *incx,cpx8 *y,int *incy,float *c,float *s); 
  void csscal_(int *n,float *a,cpx8 *x,int *incx); 
  void cswap_(int *n,cpx8 *x,int *incx,cpx8 *y,int *incy); 
  int  icamax_(int *n,cpx8 *x,int *incx); 
  int  icamin_(int *n,cpx8 *x,int *incx); 
  
  double dasum_(int *n,double *x,int *incx);
  void   daxpy_(int *n,double *alpha,double *x,int *incx,double *y,int *incy);
  void   daxpyi_(int *nz,double *a,double *x,int *indx,double *y);
  void   dcopy_(int *n,double *x,int *incx,double *y,int *incy);
  double ddot_(int *n,double *x,int *incx,double *y,int *incy);
  double ddoti_(int *nz,double *x,int *indx,double *y);
  void   dgthr_(int *nz,double *y,double *x,int *indx);
  void   dgthrz_(int *nz,double *y,double *x,int *indx);
  double dnrm2_(int *n,double *x,int *incx);
  void   drot_(int *n,double *x,int *incx,double *y,int *incy,double *c,double *s);
  void   drotg_(double *a,double *b,double *c,double *s);
  void   droti_(int *nz,double *x,int *indx,double *y,double *c,double *s);
  void   drotm_(int *n,double *x,int *incx,double *y,int *incy,double *param);
  void   drotmg_(double *d1,double *d2,double *x1,double *y1,double *param);
  void   dscal_(int *n,double *a,double *x,int *incx);
  void   dsctr_(int *nz,double *x,int *indx,double *y);
  void   dswap_(int *n,double *x,int *incx,double *y,int *incy);
  double dzasum_(int *n,cpx16 *x,int *incx); 
  double dznrm2_(int *n,cpx16 *x,int *incx); 
  int    idamax_(int *n,double *x,int *incx);
  int    idamin_(int *n,double *x,int *incx);
  
  void zaxpy_(int *n,cpx16 *alpha,cpx16 *x,int *incx,cpx16 *y,int *incy); 
  void zaxpyi_(int *nz,cpx16 *a,cpx16 *x,int *indx,cpx16 *y); 
  void zcopy_(int *n,cpx16 *x,int *incx,cpx16 *y,int *incy); 
  void zdotc_(cpx16 *pres,int *n,cpx16 *x,int *incx,cpx16 *y,int *incy); 
  void zdotci_(cpx16 *pres,int *nz,cpx16 *x,int *indx,cpx16 *y); 
  void zdotu_(cpx16 *pres,int *n,cpx16 *x,int *incx,cpx16 *y,int *incy); 
  void zdotui_(cpx16 *pres,int *nz,cpx16 *x,int *indx,cpx16 *y); 
  void zdrot_(int *n,cpx16 *x,int *incx,cpx16 *y,int *incy,double *c,double *s); 
  void zdscal_(int *n,double *a,cpx16 *x,int *incx); 
  void zgthr_(int *nz,cpx16 *y,cpx16 *x,int *indx); 
  void zgthrz_(int *nz,cpx16 *y,cpx16 *x,int *indx); 
  void zrotg_(cpx16 *a,cpx16 *b,double *c,cpx16 *s); 
  void zscal_(int *n,cpx16 *a,cpx16 *x,int *incx); 
  void zsctr_(int *nz,cpx16 *x,int *indx,cpx16 *y); 
  void zswap_(int *n,cpx16 *x,int *incx,cpx16 *y,int *incy); 
  int  izamax_(int *n,cpx16 *x,int *incx); 
  int  izamin_(int *n,cpx16 *x,int *incx); 
  
  /* blas level2 */
  
  void sgbmv_(char *trans,int *m,int *n,int *kl,int *ku,float *alpha,float *a,int *lda,float *x,int *incx,float *beta,float *y,int *incy);
  void sgemv_(char *trans,int *m,int *n,float *alpha,float *a,int *lda,float *x,int *incx,float *beta,float *y,int *incy);
  void sger_(int *m,int *n,float *alpha,float *x,int *incx,float *y,int *incy,float *a,int *lda);
  void ssbmv_(char *uplo,int *n,int *k,float *alpha,float *a,int *lda,float *x,int *incx,float *beta,float *y,int *incy);
  void sspmv_(char *uplo,int *n,float *alpha,float *ap,float *x,int *incx,float *beta,float *y,int *incy);
  void sspr_(char *uplo,int *n,float *alpha,float *x,int *incx,float *ap);
  void sspr2_(char *uplo,int *n,float *alpha,float *x,int *incx,float *y,int *incy,float *ap);
  void ssymv_(char *uplo,int *n,float *alpha,float *a,int *lda,float *x,int *incx,float *beta,float *y,int *incy);
  void ssyr_(char *uplo,int *n,float *alpha,float *x,int *incx,float *a,int *lda);
  void ssyr2_(char *uplo,int *n,float *alpha,float *x,int *incx,float *y,int *incy,float *a,int *lda);
  void stbmv_(char *uplo,char *trans,char *diag,int *n,int *k,float *a,int *lda,float *x,int *incx);
  void stbsv_(char *uplo,char *trans,char *diag,int *n,int *k,float *a,int *lda,float *x,int *incx);
  void stpmv_(char *uplo,char *trans,char *diag,int *n,float *ap,float *x,int *incx);
  void stpsv_(char *uplo,char *trans,char *diag,int *n,float *ap,float *x,int *incx);
  void strmv_(char *uplo,char *transa,char *diag,int *n,float *a,int *lda,float *b,int *incx);
  void strsv_(char *uplo,char *trans,char *diag,int *n,float *a,int *lda,float *x,int *incx);
  
  void cgbmv_(char *trans,int *m,int *n,int *kl,int *ku,cpx8 *alpha,cpx8 *a,int *lda,cpx8 *x,int *incx,cpx8 *beta,cpx8 *y,int *incy); 
  void cgemv_(char *trans,int *m,int *n,cpx8 *alpha,cpx8 *a,int *lda,cpx8 *x,int *incx,cpx8 *beta,cpx8 *y,int *incy); 
  void cgerc_(int *m,int *n,cpx8 *alpha,cpx8 *x,int *incx,cpx8 *y,int *incy,cpx8 *a,int *lda); 
  void cgeru_(int *m,int *n,cpx8 *alpha,cpx8 *x,int *incx,cpx8 *y,int *incy,cpx8 *a,int *lda); 
  void chbmv_(char *uplo,int *n,int *k,cpx8 *alpha,cpx8 *a,int *lda,cpx8 *x,int *incx,cpx8 *beta,cpx8 *y,int *incy); 
  void chemv_(char *uplo,int *n,cpx8 *alpha,cpx8 *a,int *lda,cpx8 *x,int *incx,cpx8 *beta,cpx8 *y,int *incy); 
  void cher_(char *uplo,int *n,float *alpha,cpx8 *x,int *incx,cpx8 *a,int *lda); 
  void cher2_(char *uplo,int *n,cpx8 *alpha,cpx8 *x,int *incx,cpx8 *y,int *incy,cpx8 *a,int *lda); 
  void chpmv_(char *uplo,int *n,cpx8 *alpha,cpx8 *ap,cpx8 *x,int *incx,cpx8 *beta,cpx8 *y,int *incy); 
  void chpr_(char *uplo,int *n,float *alpha,cpx8 *x,int *incx,cpx8 *ap); 
  void chpr2_(char *uplo,int *n,cpx8 *alpha,cpx8 *x,int *incx,cpx8 *y,int *incy,cpx8 *ap); 
  void ctbmv_(char *uplo,char *trans,char *diag,int *n,int *k,cpx8 *a,int *lda,cpx8 *x,int *incx); 
  void ctbsv_(char *uplo,char *trans,char *diag,int *n,int *k,cpx8 *a,int *lda,cpx8 *x,int *incx); 
  void ctpmv_(char *uplo,char *trans,char *diag,int *n,cpx8 *ap,cpx8 *x,int *incx); 
  void ctpsv_(char *uplo,char *trans,char *diag,int *n,cpx8 *ap,cpx8 *x,int *incx); 
  void ctrmv_(char *uplo,char *transa,char *diag,int *n,cpx8 *a,int *lda,cpx8 *b,int *incx); 
  void ctrsv_(char *uplo,char *trans,char *diag,int *n,cpx8 *a,int *lda,cpx8 *x,int *incx); 
  
  void dgbmv_(char *trans,int *m,int *n,int *kl,int *ku,double *alpha,double *a,int *lda,double *x,int *incx,double *beta,double *y,int *incy);
  void dgemv_(char *trans,int *m,int *n,double *alpha,double *a,int *lda,double *x,int *incx,double *beta,double *y,int *incy);
  void dger_(int *m,int *n,double *alpha,double *x,int *incx,double *y,int *incy,double *a,int *lda);
  void dsbmv_(char *uplo,int *n,int *k,double *alpha,double *a,int *lda,double *x,int *incx,double *beta,double *y,int *incy);
  void dspmv_(char *uplo,int *n,double *alpha,double *ap,double *x,int *incx,double *beta,double *y,int *incy);
  void dspr_(char *uplo,int *n,double *alpha,double *x,int *incx,double *ap);
  void dspr2_(char *uplo,int *n,double *alpha,double *x,int *incx,double *y,int *incy,double *ap);
  void dsymv_(char *uplo,int *n,double *alpha,double *a,int *lda,double *x,int *incx,double *beta,double *y,int *incy);
  void dsyr_(char *uplo,int *n,double *alpha,double *x,int *incx,double *a,int *lda);
  void dsyr2_(char *uplo,int *n,double *alpha,double *x,int *incx,double *y,int *incy,double *a,int *lda);
  void dtbmv_(char *uplo,char *trans,char *diag,int *n,int *k,double *a,int *lda,double *x,int *incx);
  void dtbsv_(char *uplo,char *trans,char *diag,int *n,int *k,double *a,int *lda,double *x,int *incx);
  void dtpmv_(char *uplo,char *trans,char *diag,int *n,double *ap,double *x,int *incx);
  void dtpsv_(char *uplo,char *trans,char *diag,int *n,double *ap,double *x,int *incx);
  void dtrmv_(char *uplo,char *transa,char *diag,int *n,double *a,int *lda,double *b,int *incx);
  void dtrsv_(char *uplo,char *trans,char *diag,int *n,double *a,int *lda,double *x,int *incx);
  
  void zgbmv_(char *trans,int *m,int *n,int *kl,int *ku,cpx16 *alpha,cpx16 *a,int *lda,cpx16 *x,int *incx,cpx16 *beta,cpx16 *y,int *incy); 
  void zgemv_(char *trans,int *m,int *n,cpx16 *alpha,cpx16 *a,int *lda,cpx16 *x,int *incx,cpx16 *beta,cpx16 *y,int *incy); 
  void zgerc_(int *m,int *n,cpx16 *alpha,cpx16 *x,int *incx,cpx16 *y,int *incy,cpx16 *a,int *lda); 
  void zgeru_(int *m,int *n,cpx16 *alpha,cpx16 *x,int *incx,cpx16 *y,int *incy,cpx16 *a,int *lda); 
  void zhbmv_(char *uplo,int *n,int *k,cpx16 *alpha,cpx16 *a,int *lda,cpx16 *x,int *incx,cpx16 *beta,cpx16 *y,int *incy); 
  void zhemv_(char *uplo,int *n,cpx16 *alpha,cpx16 *a,int *lda,cpx16 *x,int *incx,cpx16 *beta,cpx16 *y,int *incy); 
  void zher_(char *uplo,int *n,double *alpha,cpx16 *x,int *incx,cpx16 *a,int *lda); 
  void zher2_(char *uplo,int *n,cpx16 *alpha,cpx16 *x,int *incx,cpx16 *y,int *incy,cpx16 *a,int *lda); 
  void zhpmv_(char *uplo,int *n,cpx16 *alpha,cpx16 *ap,cpx16 *x,int *incx,cpx16 *beta,cpx16 *y,int *incy); 
  void zhpr_(char *uplo,int *n,double *alpha,cpx16 *x,int *incx,cpx16 *ap); 
  void zhpr2_(char *uplo,int *n,cpx16 *alpha,cpx16 *x,int *incx,cpx16 *y,int *incy,cpx16 *ap); 
  void ztbmv_(char *uplo,char *trans,char *diag,int *n,int *k,cpx16 *a,int *lda,cpx16 *x,int *incx); 
  void ztbsv_(char *uplo,char *trans,char *diag,int *n,int *k,cpx16 *a,int *lda,cpx16 *x,int *incx); 
  void ztpmv_(char *uplo,char *trans,char *diag,int *n,cpx16 *ap,cpx16 *x,int *incx); 
  void ztpsv_(char *uplo,char *trans,char *diag,int *n,cpx16 *ap,cpx16 *x,int *incx); 
  void ztrmv_(char *uplo,char *transa,char *diag,int *n,cpx16 *a,int *lda,cpx16 *b,int *incx); 
  void ztrsv_(char *uplo,char *trans,char *diag,int *n,cpx16 *a,int *lda,cpx16 *x,int *incx); 
  
  /* blas level3 */
  
  void sgemm_(char *transa,char *transb,int *m,int *n,int *k,float *alpha,float *a,int *lda,float *b,int *ldb,float *beta,float *c,int *ldc);
  void ssymm_(char *side,char *uplo,int *m,int *n,float *alpha,float *a,int *lda,float *b,int *ldb,float *beta,float *c,int *ldc);
  void ssyr2k_(char *uplo,char *trans,int *n,int *k,float *alpha,float *a,int *lda,float *b,int *ldb,float *beta,float *c,int *ldc);
  void ssyrk_(char *uplo,char *trans,int *n,int *k,float *alpha,float *a,int *lda,float *beta,float *c,int *ldc);
  void strmm_(char *side,char *uplo,char *transa,char *diag,int *m,int *n,float *alpha,float *a,int *lda,float *b,int *ldb);
  void strsm_(char *side,char *uplo,char *transa,char *diag,int *m,int *n,float *alpha,float *a,int *lda,float *b,int *ldb);
  
  void cgemm_(char *transa,char *transb,int *m,int *n,int *k,cpx8 *alpha,cpx8 *a,int *lda,cpx8 *b,int *ldb,cpx8 *beta,cpx8 *c,int *ldc); 
  void chemm_(char *side,char *uplo,int *m,int *n,cpx8 *alpha,cpx8 *a,int *lda,cpx8 *b,int *ldb,cpx8 *beta,cpx8 *c,int *ldc); 
  void cher2k_(char *uplo,char *trans,int *n,int *k,cpx8 *alpha,cpx8 *a,int *lda,cpx8 *b,int *ldb,float *beta,cpx8 *c,int *ldc); 
  void cherk_(char *uplo,char *trans,int *n,int *k,float *alpha,cpx8 *a,int *lda,float *beta,cpx8 *c,int *ldc); 
  void csymm_(char *side,char *uplo,int *m,int *n,cpx8 *alpha,cpx8 *a,int *lda,cpx8 *b,int *ldb,cpx8 *beta,cpx8 *c,int *ldc); 
  void csyr2k_(char *uplo,char *trans,int *n,int *k,cpx8 *alpha,cpx8 *a,int *lda,cpx8 *b,int *ldb,cpx8 *beta,cpx8 *c,int *ldc); 
  void csyrk_(char *uplo,char *trans,int *n,int *k,cpx8 *alpha,cpx8 *a,int *lda,cpx8 *beta,cpx8 *c,int *ldc); 
  void ctrmm_(char *side,char *uplo,char *transa,char *diag,int *m,int *n,cpx8 *alpha,cpx8 *a,int *lda,cpx8 *b,int *ldb); 
  void ctrsm_(char *side,char *uplo,char *transa,char *diag,int *m,int *n,cpx8 *alpha,cpx8 *a,int *lda,cpx8 *b,int *ldb); 
  
  void dgemm_(char *transa,char *transb,int *m,int *n,int *k,double *alpha,double *a,int *lda,double *b,int *ldb,double *beta,double *c,int *ldc);
  void dsymm_(char *side,char *uplo,int *m,int *n,double *alpha,double *a,int *lda,double *b,int *ldb,double *beta,double *c,int *ldc);
  void dsyr2k_(char *uplo,char *trans,int *n,int *k,double *alpha,double *a,int *lda,double *b,int *ldb,double *beta,double *c,int *ldc);
  void dsyrk_(char *uplo,char *trans,int *n,int *k,double *alpha,double *a,int *lda,double *beta,double *c,int *ldc);
  void dtrmm_(char *side,char *uplo,char *transa,char *diag,int *m,int *n,double *alpha,double *a,int *lda,double *b,int *ldb);
  void dtrsm_(char *side,char *uplo,char *transa,char *diag,int *m,int *n,double *alpha,double *a,int *lda,double *b,int *ldb);
  
  void zgemm_(char *transa,char *transb,int *m,int *n,int *k,cpx16 *alpha,cpx16 *a,int *lda,cpx16 *b,int *ldb,cpx16 *beta,cpx16 *c,int *ldc); 
  void zgemm_(char *transa,char *transb,int *m,int *n,int *k,cpx16 *alpha,cpx16 *a,int *lda,cpx16 *b,int *ldb,cpx16 *beta,cpx16 *c,int *ldc); 
  void zhemm_(char *side,char *uplo,int *m,int *n,cpx16 *alpha,cpx16 *a,int *lda,cpx16 *b,int *ldb,cpx16 *beta,cpx16 *c,int *ldc); 
  void zher2k_(char *uplo,char *trans,int *n,int *k,cpx16 *alpha,cpx16 *a,int *lda,cpx16 *b,int *ldb,double *beta,cpx16 *c,int *ldc); 
  void zherk_(char *uplo,char *trans,int *n,int *k,double *alpha,cpx16 *a,int *lda,double *beta,cpx16 *c,int *ldc); 
  void zsymm_(char *side,char *uplo,int *m,int *n,cpx16 *alpha,cpx16 *a,int *lda,cpx16 *b,int *ldb,cpx16 *beta,cpx16 *c,int *ldc); 
  void zsyr2k_(char *uplo,char *trans,int *n,int *k,cpx16 *alpha,cpx16 *a,int *lda,cpx16 *b,int *ldb,cpx16 *beta,cpx16 *c,int *ldc); 
  void zsyrk_(char *uplo,char *trans,int *n,int *k,cpx16 *alpha,cpx16 *a,int *lda,cpx16 *beta,cpx16 *c,int *ldc); 
  void ztrmm_(char *side,char *uplo,char *transa,char *diag,int *m,int *n,cpx16 *alpha,cpx16 *a,int *lda,cpx16 *b,int *ldb); 
  void ztrsm_(char *side,char *uplo,char *transa,char *diag,int *m,int *n,cpx16 *alpha,cpx16 *a,int *lda,cpx16 *b,int *ldb); 
}

#endif

