#ifndef _sf_fblas_h_
#define _sf_fblas_h_

#include "komplex.h"

#ifndef cpx8
#define cpx8 sf_complex
#endif

#ifndef cpx16
#define cpx16 sf_double_complex
#endif
  
double sasum_(int *n,float *x,int *incx);
int  saxpy_(int *n,float *alpha,float *x,int *incx,float *y,int *incy);
void  saxpyi_(int *nz,float *a,float *x,int *indx,float *y);
double scasum_(int *n,void *x,int *incx); 
double scnrm2_(int *n,void *x,int *incx); 
int  scopy_(int *n,float *x,int *incx,float *y,int *incy);
double sdot_(int *n,float *x,int *incx,float *y,int *incy);
float sdoti_(int *nz,float *x,int *indx,float *y);
void  sgthr_(int *nz,float *y,float *x,int *indx);
void  sgthrz_(int *nz,float *y,float *x,int *indx);
double snrm2_(int *n,float *x,int *incx);
int  srot_(int *n,float *x,int *incx,float *y,int *incy,float *c,float *s);
int  srotg_(float *a,float *b,float *c,float *s);
void  sroti_(int *nz,float *x,int *indx,float *y,float *c,float *s);
int  srotm_(int *n,float *x,int *incx,float *y,int *incy,float *param);
int  srotmg_(float *d1,float *d2,float *x1,float *y1,float *param);
int  sscal_(int *n,float *a,float *x,int *incx);
void  ssctr_(int *nz,float *x,int *indx,float *y);
int  sswap_(int *n,float *x,int *incx,float *y,int *incy);
int   isamax_(int *n,float *x,int *incx);
int   isamin_(int *n,float *x,int *incx);
  
int caxpy_(int *n,void *alpha,void *x,int *incx,void *y,int *incy); 
void caxpyi_(int *nz,cpx8 *a,cpx8 *x,int *indx,cpx8 *y); 
int ccopy_(int *n,void *x,int *incx,void *y,int *incy); 
void cdotc_(void *pres,int *n,void *x,int *incx,void *y,int *incy); 
void cdotci_(cpx8 *pres,int *nz,cpx8 *x,int *indx,cpx8 *y); 
void cdotu_(void *pres,int *n,void *x,int *incx,void *y,int *incy); 
void cdotui_(cpx8 *pres,int *nz,cpx8 *x,int *indx,cpx8 *y); 
void cgthr_(int *nz,cpx8 *y,cpx8 *x,int *indx); 
void cgthrz_(int *nz,cpx8 *y,cpx8 *x,int *indx); 
int crotg_(void *a,void *b,float *c,void *s); 
int cscal_(int *n,void *a,void *x,int *incx); 
void csctr_(int *nz,cpx8 *x,int *indx,cpx8 *y); 
int csrot_(int *n,void *x,int *incx,void *y,int *incy,float *c,float *s); 
int csscal_(int *n,float *a,void *x,int *incx); 
int cswap_(int *n,void *x,int *incx,void *y,int *incy); 
int  icamax_(int *n,void *x,int *incx); 
int  icamin_(int *n,cpx8 *x,int *incx); 
  
double dasum_(int *n,double *x,int *incx);
int   daxpy_(int *n,double *alpha,double *x,int *incx,double *y,int *incy);
void   daxpyi_(int *nz,double *a,double *x,int *indx,double *y);
int   dcopy_(int *n,double *x,int *incx,double *y,int *incy);
double ddot_(int *n,double *x,int *incx,double *y,int *incy);
double ddoti_(int *nz,double *x,int *indx,double *y);
void   dgthr_(int *nz,double *y,double *x,int *indx);
void   dgthrz_(int *nz,double *y,double *x,int *indx);
double dnrm2_(int *n,double *x,int *incx);
int   drot_(int *n,double *x,int *incx,double *y,int *incy,double *c,double *s);
int   drotg_(double *a,double *b,double *c,double *s);
void   droti_(int *nz,double *x,int *indx,double *y,double *c,double *s);
int   drotm_(int *n,double *x,int *incx,double *y,int *incy,double *param);
int   drotmg_(double *d1,double *d2,double *x1,double *y1,double *param);
int   dscal_(int *n,double *a,double *x,int *incx);
void   dsctr_(int *nz,double *x,int *indx,double *y);
int   dswap_(int *n,double *x,int *incx,double *y,int *incy);
double dzasum_(int *n,void *x,int *incx); 
double dznrm2_(int *n,void *x,int *incx); 
int    idamax_(int *n,double *x,int *incx);
int    idamin_(int *n,double *x,int *incx);
  
int zaxpy_(int *n,void *alpha,void *x,int *incx,void *y,int *incy); 
void zaxpyi_(int *nz,cpx16 *a,cpx16 *x,int *indx,cpx16 *y); 
int zcopy_(int *n,void *x,int *incx,void *y,int *incy); 
void zdotc_(void *pres,int *n,void *x,int *incx,void *y,int *incy); 
void zdotci_(cpx16 *pres,int *nz,cpx16 *x,int *indx,cpx16 *y); 
void zdotu_(void *pres,int *n,void *x,int *incx,void *y,int *incy); 
void zdotui_(cpx16 *pres,int *nz,cpx16 *x,int *indx,cpx16 *y); 
int zdrot_(int *n,void *x,int *incx,void *y,int *incy,double *c,double *s); 
int zdscal_(int *n,double *a,void *x,int *incx); 
void zgthr_(int *nz,cpx16 *y,cpx16 *x,int *indx); 
void zgthrz_(int *nz,cpx16 *y,cpx16 *x,int *indx); 
int zrotg_(void *a,void *b,double *c,void *s); 
int zscal_(int *n,void *a,void *x,int *incx); 
void zsctr_(int *nz,cpx16 *x,int *indx,cpx16 *y); 
int zswap_(int *n,void *x,int *incx,void *y,int *incy); 
int  izamax_(int *n,void *x,int *incx); 
int  izamin_(int *n,cpx16 *x,int *incx); 
  
/* blas level2 */
  
int sgbmv_(char *trans,int *m,int *n,int *kl,int *ku,float *alpha,float *a,int *lda,float *x,int *incx,float *beta,float *y,int *incy);
int sgemv_(char *trans,int *m,int *n,float *alpha,float *a,int *lda,float *x,int *incx,float *beta,float *y,int *incy);
int sger_(int *m,int *n,float *alpha,float *x,int *incx,float *y,int *incy,float *a,int *lda);
int ssbmv_(char *uplo,int *n,int *k,float *alpha,float *a,int *lda,float *x,int *incx,float *beta,float *y,int *incy);
int sspmv_(char *uplo,int *n,float *alpha,float *ap,float *x,int *incx,float *beta,float *y,int *incy);
int sspr_(char *uplo,int *n,float *alpha,float *x,int *incx,float *ap);
int sspr2_(char *uplo,int *n,float *alpha,float *x,int *incx,float *y,int *incy,float *ap);
int ssymv_(char *uplo,int *n,float *alpha,float *a,int *lda,float *x,int *incx,float *beta,float *y,int *incy);
int ssyr_(char *uplo,int *n,float *alpha,float *x,int *incx,float *a,int *lda);
int ssyr2_(char *uplo,int *n,float *alpha,float *x,int *incx,float *y,int *incy,float *a,int *lda);
int stbmv_(char *uplo,char *trans,char *diag,int *n,int *k,float *a,int *lda,float *x,int *incx);
int stbsv_(char *uplo,char *trans,char *diag,int *n,int *k,float *a,int *lda,float *x,int *incx);
int stpmv_(char *uplo,char *trans,char *diag,int *n,float *ap,float *x,int *incx);
int stpsv_(char *uplo,char *trans,char *diag,int *n,float *ap,float *x,int *incx);
int strmv_(char *uplo,char *transa,char *diag,int *n,float *a,int *lda,float *b,int *incx);
int strsv_(char *uplo,char *trans,char *diag,int *n,float *a,int *lda,float *x,int *incx);
  
int cgbmv_(char *trans,int *m,int *n,int *kl,int *ku,void *alpha,void *a,int *lda,void *x,int *incx,void *beta,void *y,int *incy); 
int cgemv_(char *trans,int *m,int *n,void *alpha,void *a,int *lda,void *x,int *incx,void *beta,void *y,int *incy); 
int cgerc_(int *m,int *n,void *alpha,void *x,int *incx,void *y,int *incy,void *a,int *lda); 
int cgeru_(int *m,int *n,void *alpha,void *x,int *incx,void *y,int *incy,void *a,int *lda); 
int chbmv_(char *uplo,int *n,int *k,void *alpha,void *a,int *lda,void *x,int *incx,void *beta,void *y,int *incy); 
int chemv_(char *uplo,int *n,void *alpha,void *a,int *lda,void *x,int *incx,void *beta,void *y,int *incy); 
int cher_(char *uplo,int *n,float *alpha,void *x,int *incx,void *a,int *lda); 
int cher2_(char *uplo,int *n,void *alpha,void *x,int *incx,void *y,int *incy,void *a,int *lda); 
int chpmv_(char *uplo,int *n,void *alpha,void *ap,void *x,int *incx,void *beta,void *y,int *incy); 
int chpr_(char *uplo,int *n,float *alpha,void *x,int *incx,void *ap); 
int chpr2_(char *uplo,int *n,void *alpha,void *x,int *incx,void *y,int *incy,void *ap); 
int ctbmv_(char *uplo,char *trans,char *diag,int *n,int *k,void *a,int *lda,void *x,int *incx); 
int ctbsv_(char *uplo,char *trans,char *diag,int *n,int *k,void *a,int *lda,void *x,int *incx); 
int ctpmv_(char *uplo,char *trans,char *diag,int *n,void *ap,void *x,int *incx); 
int ctpsv_(char *uplo,char *trans,char *diag,int *n,void *ap,void *x,int *incx); 
int ctrmv_(char *uplo,char *transa,char *diag,int *n,void *a,int *lda,void *b,int *incx); 
int ctrsv_(char *uplo,char *trans,char *diag,int *n,void *a,int *lda,void *x,int *incx); 
  
int dgbmv_(char *trans,int *m,int *n,int *kl,int *ku,double *alpha,double *a,int *lda,double *x,int *incx,double *beta,double *y,int *incy);
int dgemv_(char *trans,int *m,int *n,double *alpha,double *a,int *lda,double *x,int *incx,double *beta,double *y,int *incy);
int dger_(int *m,int *n,double *alpha,double *x,int *incx,double *y,int *incy,double *a,int *lda);
int dsbmv_(char *uplo,int *n,int *k,double *alpha,double *a,int *lda,double *x,int *incx,double *beta,double *y,int *incy);
int dspmv_(char *uplo,int *n,double *alpha,double *ap,double *x,int *incx,double *beta,double *y,int *incy);
int dspr_(char *uplo,int *n,double *alpha,double *x,int *incx,double *ap);
int dspr2_(char *uplo,int *n,double *alpha,double *x,int *incx,double *y,int *incy,double *ap);
int dsymv_(char *uplo,int *n,double *alpha,double *a,int *lda,double *x,int *incx,double *beta,double *y,int *incy);
int dsyr_(char *uplo,int *n,double *alpha,double *x,int *incx,double *a,int *lda);
int dsyr2_(char *uplo,int *n,double *alpha,double *x,int *incx,double *y,int *incy,double *a,int *lda);
int dtbmv_(char *uplo,char *trans,char *diag,int *n,int *k,double *a,int *lda,double *x,int *incx);
int dtbsv_(char *uplo,char *trans,char *diag,int *n,int *k,double *a,int *lda,double *x,int *incx);
int dtpmv_(char *uplo,char *trans,char *diag,int *n,double *ap,double *x,int *incx);
int dtpsv_(char *uplo,char *trans,char *diag,int *n,double *ap,double *x,int *incx);
int dtrmv_(char *uplo,char *transa,char *diag,int *n,double *a,int *lda,double *b,int *incx);
int dtrsv_(char *uplo,char *trans,char *diag,int *n,double *a,int *lda,double *x,int *incx);
  
int zgbmv_(char *trans,int *m,int *n,int *kl,int *ku,void *alpha,void *a,int *lda,void *x,int *incx,void *beta,void *y,int *incy); 
int zgemv_(char *trans,int *m,int *n,void *alpha,void *a,int *lda,void *x,int *incx,void *beta,void *y,int *incy); 
int zgerc_(int *m,int *n,void *alpha,void *x,int *incx,void *y,int *incy,void *a,int *lda); 
int zgeru_(int *m,int *n,void *alpha,void *x,int *incx,void *y,int *incy,void *a,int *lda); 
int zhbmv_(char *uplo,int *n,int *k,void *alpha,void *a,int *lda,void *x,int *incx,void *beta,void *y,int *incy); 
int zhemv_(char *uplo,int *n,void *alpha,void *a,int *lda,void *x,int *incx,void *beta,void *y,int *incy); 
int zher_(char *uplo,int *n,double *alpha,void *x,int *incx,void *a,int *lda); 
int zher2_(char *uplo,int *n,void *alpha,void *x,int *incx,void *y,int *incy,void *a,int *lda); 
int zhpmv_(char *uplo,int *n,void *alpha,void *ap,void *x,int *incx,void *beta,void *y,int *incy); 
int zhpr_(char *uplo,int *n,double *alpha,void *x,int *incx,void *ap); 
int zhpr2_(char *uplo,int *n,void *alpha,void *x,int *incx,void *y,int *incy,void *ap); 
int ztbmv_(char *uplo,char *trans,char *diag,int *n,int *k,void *a,int *lda,void *x,int *incx); 
int ztbsv_(char *uplo,char *trans,char *diag,int *n,int *k,void *a,int *lda,void *x,int *incx); 
int ztpmv_(char *uplo,char *trans,char *diag,int *n,void *ap,void *x,int *incx); 
int ztpsv_(char *uplo,char *trans,char *diag,int *n,void *ap,void *x,int *incx); 
int ztrmv_(char *uplo,char *transa,char *diag,int *n,void *a,int *lda,void *b,int *incx); 
int ztrsv_(char *uplo,char *trans,char *diag,int *n,void *a,int *lda,void *x,int *incx); 
  
/* blas level3 */
  
int sgemm_(char *transa,char *transb,int *m,int *n,int *k,float *alpha,float *a,int *lda,float *b,int *ldb,float *beta,float *c,int *ldc);
int ssymm_(char *side,char *uplo,int *m,int *n,float *alpha,float *a,int *lda,float *b,int *ldb,float *beta,float *c,int *ldc);
int ssyr2k_(char *uplo,char *trans,int *n,int *k,float *alpha,float *a,int *lda,float *b,int *ldb,float *beta,float *c,int *ldc);
int ssyrk_(char *uplo,char *trans,int *n,int *k,float *alpha,float *a,int *lda,float *beta,float *c,int *ldc);
int strmm_(char *side,char *uplo,char *transa,char *diag,int *m,int *n,float *alpha,float *a,int *lda,float *b,int *ldb);
int strsm_(char *side,char *uplo,char *transa,char *diag,int *m,int *n,float *alpha,float *a,int *lda,float *b,int *ldb);
  
int cgemm_(char *transa,char *transb,int *m,int *n,int *k,void *alpha,void *a,int *lda,void *b,int *ldb,void *beta,void *c,int *ldc); 
int chemm_(char *side,char *uplo,int *m,int *n,void *alpha,void *a,int *lda,void *b,int *ldb,void *beta,void *c,int *ldc); 
int cher2k_(char *uplo,char *trans,int *n,int *k,void *alpha,void *a,int *lda,void *b,int *ldb,float *beta,void *c,int *ldc); 
int cherk_(char *uplo,char *trans,int *n,int *k,float *alpha,void *a,int *lda,float *beta,void *c,int *ldc); 
int csymm_(char *side,char *uplo,int *m,int *n,void *alpha,void *a,int *lda,void *b,int *ldb,void *beta,void *c,int *ldc); 
int csyr2k_(char *uplo,char *trans,int *n,int *k,void *alpha,void *a,int *lda,void *b,int *ldb,void *beta,void *c,int *ldc); 
int csyrk_(char *uplo,char *trans,int *n,int *k,void *alpha,void *a,int *lda,void *beta,void *c,int *ldc); 
int ctrmm_(char *side,char *uplo,char *transa,char *diag,int *m,int *n,void *alpha,void *a,int *lda,void *b,int *ldb); 
int ctrsm_(char *side,char *uplo,char *transa,char *diag,int *m,int *n,void *alpha,void *a,int *lda,void *b,int *ldb); 
  
int dgemm_(char *transa,char *transb,int *m,int *n,int *k,double *alpha,double *a,int *lda,double *b,int *ldb,double *beta,double *c,int *ldc);
int dsymm_(char *side,char *uplo,int *m,int *n,double *alpha,double *a,int *lda,double *b,int *ldb,double *beta,double *c,int *ldc);
int dsyr2k_(char *uplo,char *trans,int *n,int *k,double *alpha,double *a,int *lda,double *b,int *ldb,double *beta,double *c,int *ldc);
int dsyrk_(char *uplo,char *trans,int *n,int *k,double *alpha,double *a,int *lda,double *beta,double *c,int *ldc);
int dtrmm_(char *side,char *uplo,char *transa,char *diag,int *m,int *n,double *alpha,double *a,int *lda,double *b,int *ldb);
int dtrsm_(char *side,char *uplo,char *transa,char *diag,int *m,int *n,double *alpha,double *a,int *lda,double *b,int *ldb);
  
int zgemm_(char *transa,char *transb,int *m,int *n,int *k,void *alpha,void *a,int *lda,void *b,int *ldb,void *beta,void *c,int *ldc); 
int zgemm_(char *transa,char *transb,int *m,int *n,int *k,void *alpha,void *a,int *lda,void *b,int *ldb,void *beta,void *c,int *ldc); 
int zhemm_(char *side,char *uplo,int *m,int *n,void *alpha,void *a,int *lda,void *b,int *ldb,void *beta,void *c,int *ldc); 
int zher2k_(char *uplo,char *trans,int *n,int *k,void *alpha,void *a,int *lda,void *b,int *ldb,double *beta,void *c,int *ldc); 
int zherk_(char *uplo,char *trans,int *n,int *k,double *alpha,void *a,int *lda,double *beta,void *c,int *ldc); 
int zsymm_(char *side,char *uplo,int *m,int *n,void *alpha,void *a,int *lda,void *b,int *ldb,void *beta,void *c,int *ldc); 
int zsyr2k_(char *uplo,char *trans,int *n,int *k,void *alpha,void *a,int *lda,void *b,int *ldb,void *beta,void *c,int *ldc); 
int zsyrk_(char *uplo,char *trans,int *n,int *k,void *alpha,void *a,int *lda,void *beta,void *c,int *ldc); 
int ztrmm_(char *side,char *uplo,char *transa,char *diag,int *m,int *n,void *alpha,void *a,int *lda,void *b,int *ldb); 
int ztrsm_(char *side,char *uplo,char *transa,char *diag,int *m,int *n,void *alpha,void *a,int *lda,void *b,int *ldb); 


#endif

