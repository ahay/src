#ifndef _UMIN_TR_MAT
#define _UMIN_TR_MAT

#include "rnop.hh"
#include "alg.hh"

extern "C" {
  typedef struct { float r, i; } f2c_complex; 
  typedef struct { double r, i; } f2c_dblcomplex; 

  int sgelsy_(long *m, long *n, long *nrhs, float *a, long * lda, float * b, long * ldb, long * jpvt, float * rcond, long * rank, float * work, long * lwork, long * info);
  int dgelsy_(long *m, long *n, long *nrhs, double *a, long * lda, double * b, long * ldb, long * jpvt, double * rcond, long * rank, double * work, long * lwork, long * info);
  int cgelsy_(long *m, long *n, long *nrhs, f2c_complex *a, long * lda, f2c_complex * b, long * ldb, long * jpvt, float * rcond, long * rank, f2c_complex * work, long * lwork, float * rwork, long * info);
  int zgelsy_(long *m, long *n, long *nrhs, f2c_dblcomplex *a, long * lda, f2c_dblcomplex * b, long * ldb, long * jpvt, double * rcond, long * rank, f2c_dblcomplex * work, long * lwork, double * rwork, long * info);
}

/*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK)) */
/*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK. */

/*  LWORK   (input) INTEGER */
/*          The dimension of the array WORK. */
/*          The unblocked strategy requires that: */
/*             LWORK >= MAX( MN+3*N+1, 2*MN+NRHS ), */
/*          where MN = min( M, N ). */
/*          The block algorithm requires that: */
/*             LWORK >= MAX( MN+2*N+NB*(N+1), 2*MN+NB*NRHS ), */
/*          where NB is an upper bound on the blocksize returned */
/*          by ILAENV for the routines SGEQP3, STZRZF, STZRQF, SORMQR, */
/*          and SORMRZ. */

// assumptions: (1) expect problems to be well-conditioned, so set rcond=1
// (3) set up for one rhs
// (2) assume that block strategy is used and N \le M, set NB=N thus
// LWORK = MAX(MN+3*N+N^2, 2*MN+N+1), WORK=LWORK 

namespace RVLUmin {

  using namespace RVL;
  using namespace RVLAlg;

  /** LAPACK-derived least squares solver. Uses the xgelsy family of
      algorithms. Since the various numeric types are opaque to C++,
      must provide four separate template specializations.

      This version is set up to conform to the CLAPACK calling
      conventions, and must be linked with CLAPACK. It is probably
      necessary to modify the interfaces slightly for Fortran-based
      vendor or compiled LAPACKs.

      assumptions: 

      <ul>
      <li>expect problems to be well-conditioned, so set rcond=1</li>
      <li>set up for one rhs</li>
      <li>assume that block strategy is used and N <= M, set NB=N thus
      LWORK = MAX(MN+3*N+N^2, 2*MN+N+1), WORK=LWORK </li>
      </ul>

      notes:

      <ul> <li>CLAPACK integer type is long - this is same as int
      under linux, but NOT under OS-X, so use long as the integer type
      for all arguments and temporaries</li> 

      <li>Since LAPACK overwrites the storage for both the matrix
      (input matrix on call, QR factors on return) and the solution
      vector (right-hand side on call, solution on return), these must
      be copied. This function collection makes a temporary for each -
      for the solution, long enough to hold either the RHS or the
      solution - then copies as necessary. This waste of storage would
      have to occur either here or in the driver, and lets the input
      matrix be passed as consts </li> 

      <li>Even if you didn't wanted to keep the LAPACK calling
      convention (overwriting inputs), you couldn't for the complex
      case, because the CLAPACK/F2C complex types are incompatible
      with std::complex.  </li>

      </ul>
  */
  template<typename T, typename A> 
  void LAPACK_LS(T * xdata,
		 T const * Adata,
		 T const * bdata,
		 A * rcond, long * rank, long m, long n, long lwork, long * info);
  
  template<>
  void LAPACK_LS<float,float>(float * xdata,
			      float const * Adata,
			      float const * bdata,
			      float * rcond, long * rank, long m, long n, long lwork, long * info) {
		    
    float * work = new float[lwork];
    for (int i=0;i<lwork;i++) work[i]=0.0f;
    *info = 0;
    long * jpvt = new long[n];
    for (int i=0;i<n;i++) jpvt[i]=0;
    float * tmp = new float[m];
    for (int i=0;i<m;i++) tmp[i]=bdata[i];
    float * tmpa = new float[m*n];
    for (int j=0;j<n;j++) {
      for (int i=0;i<m;i++) {
	tmpa[i+j*m]=Adata[i+j*m];
      }
    }
    long nrhs = 1;
    long lda = m;
    long ldb = m;
    sgelsy_(&m, &n, &nrhs, tmpa, &lda, tmp, &ldb, jpvt, rcond, rank, work, &lwork, info);

    if (*info==0) {
      for (int i=0;i<n;i++) {
	xdata[i]=tmp[i];
      }
    }

    delete [] jpvt;
    delete [] work;
    delete [] tmpa;
    delete [] tmp;
  }

  template<>
  void LAPACK_LS<double,double>(double * xdata,
			      double const * Adata,
			      double const * bdata,
			      double * rcond, long * rank, long m, long n, long lwork, long * info) {
		    
    double * work = new double[lwork];
    for (int i=0;i<lwork;i++) work[i]=0.0;
    *info = 0;
    long * jpvt = new long[n];
    for (int i=0;i<n;i++) jpvt[i]=0;
    double * tmp = new double[m];
    for (int i=0;i<m;i++) tmp[i]=bdata[i];
    double * tmpa = new double[m*n];
    for (int j=0;j<n;j++) {
      for (int i=0;i<m;i++) {
	tmpa[i+j*m]=Adata[i+j*m];
      }
    }
    long nrhs = 1;
    long lda = m;
    long ldb = m;
    dgelsy_(&m, &n, &nrhs, tmpa, &lda, tmp, &ldb, jpvt, rcond, rank, work, &lwork, info);

    if (*info==0) {
      for (int i=0;i<n;i++) {
	xdata[i]=tmp[i];
      }
    }

    delete [] jpvt;
    delete [] work;
    delete [] tmpa;
    delete [] tmp;
  }

  template<>
  void LAPACK_LS<std::complex<float>,float>(std::complex<float> * xdata,
					    std::complex<float> const * Adata,
					    std::complex<float> const * bdata,
					    float * rcond, long * rank, long m, 
					    long n, long lwork, long * info) {
		    
    f2c_complex * work = new f2c_complex[lwork];
    for (int i=0;i<lwork;i++) { work[i].r=0.0f; work[i].i=0.0f; }
    *info = 0;
    long * jpvt = new long[n];
    for (int i=0;i<n;i++) jpvt[i]=0;
    f2c_complex * tmp = new f2c_complex[m];
    for (int i=0;i<m;i++) { tmp[i].r=real(bdata[i]); tmp[i].i=imag(bdata[i]); }
    f2c_complex * tmpa = new f2c_complex[m*n];
    float * rwork = new float[2*n];
    for (int i=0;i<2*n;i++) rwork[i]=0.0f;
    for (int j=0;j<n;j++) {
      for (int i=0;i<m;i++) {
	tmpa[i+j*m].r=real(Adata[i+j*m]);
	tmpa[i+j*m].i=imag(Adata[i+j*m]);	
      }
    }
    long nrhs = 1;
    long lda = m;
    long ldb = m;

    cgelsy_(&m, &n, &nrhs, tmpa, &lda, tmp, &ldb, jpvt, rcond, rank, work, &lwork, rwork, info);

    if (*info==0) {
      for (int j=0;j<n;j++) {
	xdata[j]=std::complex<float>(tmp[j].r,tmp[j].i);
      }
    }

    delete [] jpvt;
    delete [] work;
    delete [] tmpa;
    delete [] tmp;
    delete [] rwork;

  }

  template<>
  void LAPACK_LS<std::complex<double>,double>(std::complex<double> * xdata,
					      std::complex<double> const * Adata,
					      std::complex<double> const * bdata,
					      double * rcond, long * rank, long m, 
					      long n, long lwork, long * info) {
		    
    f2c_dblcomplex * work = new f2c_dblcomplex[lwork];
    for (int i=0;i<lwork;i++) { work[i].r=0.0; work[i].i=0.0; }
    *info = 0;
    long * jpvt = new long[n];
    for (int i=0;i<n;i++) jpvt[i]=0;
    f2c_dblcomplex * tmp = new f2c_dblcomplex[m];
    for (int i=0;i<m;i++) { tmp[i].r=real(bdata[i]); tmp[i].i=imag(bdata[i]); }
    f2c_dblcomplex * tmpa = new f2c_dblcomplex[m*n];
    double * rwork = new double[2*n];
    for (int i=0;i<2*n;i++) rwork[i]=0.0f;
    for (int j=0;j<n;j++) {
      for (int i=0;i<m;i++) {
	tmpa[i+j*m].r=real(Adata[i+j*m]);
	tmpa[i+j*m].i=imag(Adata[i+j*m]);	
      }
    }
    long nrhs = 1;
    long lda = m;
    long ldb = m;
    zgelsy_(&m, &n, &nrhs, tmpa, &lda, tmp, &ldb, jpvt, rcond, rank, work, &lwork, rwork, info);
    if (*info==0) {
      for (int i=0;i<n;i++) {
	xdata[i]=std::complex<double>(tmp[i].r,tmp[i].i);
      }
    }

    delete [] jpvt;
    delete [] work;
    delete [] tmpa;
    delete [] tmp;
    delete [] rwork;
    
  }

  /** LAPACK-based least squares solver as a local FO, with 
      the matrix stored as instance data for a GenMat
  */
  template<typename T> 
  class LSMatFO: public BinaryLocalFunctionObject<T> {

    typedef typename ScalarFieldTraits<T>::AbsType AT;

  public:

    void operator()(LocalDataContainer<T> & x,
		    LocalDataContainer<T> const & b) {
      try {
	long m = A.getNRows();
	long n = A.getNCols();
	AT rcond=0.00001;
	long rank=0;
	long lwork = max((m+n+3)*n,(3*m+1)*n+1);
	long info=0;
	LAPACK_LS<T,AT>(x.getData(),
			A.getData(),
			b.getData(),
			&rcond,&rank,m,n,lwork,&info);
	if (info<0) {
	  RVLException e;
	  e<<"Error: LAPACK_LS, float case\n";
	  e<<"SGELSY returned with nozero error flag = "<<info<<"\n";
	  throw e;
	}	
      }
      catch (RVLException & e) {
	e<<"\ncalled from LSMat<T>::run\n";
	throw e;
      }
    }

    /** main constructor builds on a GenMat */
    LSMatFO(GenMat<T> const & _A): A(_A) {}

    string getName() const { string tmp = "LSMatFO"; return tmp; }

  private:
    
    GenMat<T> const & A;

  };

  /** LAPACK-based least squares solver with trust region truncation.

      Stores mutable references to 
      <ul>
      <li>solution vector - zeroed on call, contains solution estimate on return</li>
     
      <li>scalars for residual norm and normal residual norm - on
      construction, these are set to the norm of the right-hand side
      and of the right-hand side of the normal equations,
      respectively, so the constructor can be used in a driver to
      initialize these quantities. On return, residual resp. normal
      residual norms at estimated solution.</li>

      </ul>

      Applies trust region truncation to solution estimate - thus the
      returned solution may NOT be the least squares solution, but
      rather the closest point in the trust region to it!

      Note that trust region solvers are also set up to act as
      Terminators, with the query function returning true if the trust
      radius is exceeded. Dependent Gauss-Newton iterations stop
      updating their solution estimates, and update trust radius
      instead, when this happens.
  */

  template<typename T>
  class LSMatAlg: public Algorithm, public Terminator {

    typedef typename ScalarFieldTraits<T>::AbsType atype;

  public:

    bool query() { return proj; }

    void run() {
      try {
	// initialize x
	x.zero();
	// solve least squares problem
	LSMatFO<T> mfo(A);
	x.eval(mfo,b);
	// scale to trust region boundary
	atype xnorm=x.norm();
	proj=false;
	if (xnorm>Delta) {
	  proj=true;
	  T scfac = Delta/xnorm;
	  x.scale(scfac);
	  str<<"RVLUmin::LSMatAlg::run: trust region truncation applied\n";
	  str<<"  untruncated solution norm = "<<xnorm<<" trust radius = "<<Delta<<"\n";
	}
	// compute residual
	A.applyOp(x,r);
	r.linComb(-ScalarFieldTraits<T>::One(),b);
	// note that this is Ax-b, rather than b-Ax - doesn't matter, as it's internal
	rnorm=r.norm();
	// compute gradient
	A.applyAdjOp(r,g);
	// signs cancel so g is least-squares gradient
	nrnorm=g.norm();
	/*
	str<<"LSMatAlg: LAPACK-based Least Squares solver with trust-region truncation"<<endl;
	str<<"trust radius              = "<<Delta<<endl;
	str<<"norm of LS solution       = "<<xnorm<<endl;
	str<<"norm of trunc LS solution = "<<xtrnorm<<endl;
	str<<"norm of residual          = "<<rnorm<<endl;
	str<<"norm of normal residual   = "<<nrnorm<<endl;
	*/
      }
      catch (RVLException & e) {
	e<<"\ncalled from LSMatAlg::run\n";
	throw e;
      }
    }
    
    /** note that initialization of rnorm, nrnorm is mandated -
	standard behaviour of LS solvers for inclusion in TR alg.
	
	@param _x - stored mutable reference. zeroed on call, solution estimate on return
	@param _A - stored const reference. GenMat representation of matrix linear operator
	@param _b - stored const reference. right-hand side (fit data) of least squares problem
	@param _rnorm - stored mutable reference. on call, set to norm of RHS. on return, residual = value of least squares function at estimated solution
	@param _nrnorm - stored mutable reference. on call, set to norm of normal residual. on return, normal residual at least squares estimate
	@param _Delta - trust radius. solution truncated to this size if necessary.
	@param _str - verbose output unit
     */
    LSMatAlg(Vector<T> & _x, 
	     GenMat<T> const & _A, 
	     Vector<T> const & _b,
	     atype & _rnorm,
	     atype & _nrnorm,
	     atype  _Delta,
	     ostream & _str)
      : x(_x), A(_A), r(A.getRange()), g(A.getDomain()), b(_b), rnorm(_rnorm), nrnorm(_nrnorm), Delta(_Delta), proj(false), str(_str) {
      r.copy(b);
      rnorm=r.norm();
      A.applyAdjOp(r,g);
      nrnorm=g.norm();
    }

  private:
    Vector<T> & x;
    GenMat<T> const & A;
    Vector<T> const & b;
    Vector<T> r;
    Vector<T> g;
    atype & rnorm;
    atype & nrnorm;
    atype Delta;
    mutable bool proj;
    ostream & str;
  };

  /** Factory class for LSMatAlg. Intended as policy for least squares
      / trust-region solver in Steihaug-Toint algorithm
      implementation, hence takes OperatorEvaluation. The Operator is
      extracted from this on construction, and an attempt is made to
      cast the Operator to a GenOp, i.e. an Operator with Derivative
      represented by a GenMat. The resulting GenMat is passed to the
      build method to implement an LSMatAlg factory.

      Conforms to specifications described in TRGNAlg docs.
  */

  template<typename T>
  class LSMatPolicy {
    
    typedef typename ScalarFieldTraits<T>::AbsType atype;
    
  public:
    
    /** Sets up Gauss-Newton problem, verifies that the derivative is
	actually represented by a GenMat, then builds solver. 
    */
    LSMatAlg<T> * build(Vector<T> & x, 
			OperatorEvaluation<T> & opeval,
			atype & rnorm,
			atype & nrnorm,
			atype Delta,
			ostream & str) const {

      try {
	// first make sure operator is applied - want the Operator for the current point.
	Vector<T> const & rhs = opeval.getValue();
	// then get Operator 
	OpWithGenMatDeriv<T> const & gop = dynamic_cast< OpWithGenMatDeriv<T> const &> (opeval.getOp());
	GenMat<T> const & G = gop.getGenMat();
	if (G.getNRows() < G.getNCols()) {
	  RVLException e;
	  e<<"Error: LSMatPolicy::build \n";
	  e<<"GenMat object passed as 2nd arg has more cols than rows\n";
	  throw e;
	}
	return new LSMatAlg<T>(x,G,rhs,rnorm,nrnorm,Delta,str);
      }
      catch (bad_cast) {
	RVLException e;
	e<<"Error: LSMatPolicy::build \n";	
	e<<"input LinearOp not of GenMat type\n";
	throw e;
      }
      catch (RVLException & e) {
	e<<"\ncalled from LSMatPolicy::build\n";
	throw e;
      }
    }

  };
}
#endif

      
