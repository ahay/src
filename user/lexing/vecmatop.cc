#include <rsf.hh>

#include "numvec.hh"
#include "nummat.hh"
#include "numtns.hh"

#include "vecmatop.hh"

using std::cerr;
using std::min;
using std::max;
//typedef long int int;

//Y <- a M X + b Y
// ---------------------------------------------------------------------- 
int dgemm(float alpha, const FltNumMat& A, const FltNumMat& B, float beta, FltNumMat& C)
{
  assert( A.m() == C.m() );  assert( A.n() == B.m() );  assert( B.n() == C.n() );
  if(A.m()==0 || A.n()==0 || B.n()==0) { //simply scale
    for(int i=0; i<C.m(); i++)
      for(int j=0; j<C.n(); j++)
	C(i,j) = beta*C(i,j);
  } else {
    iC( dgemm(C.m(), C.n(), A.n(), alpha, A.data(), B.data(), beta, C.data()) );
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int dgemm(int m, int n, int k, float alpha, float* A, float* B, float beta, float* C)
{
  char transa = 'N';
  char transb = 'N';
  assert(m!=0 && n!=0 && k!=0);
  sgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &m, B, &k, &beta, C, &m);
  return 0;
}
//Y <- a M X + b Y
// ---------------------------------------------------------------------- 
int dgemv(float alpha, const FltNumMat& A, const FltNumVec& X, float beta, FltNumVec& Y)
{
  assert(Y.m() == A.m());
  assert(A.n() == X.m());
  if(A.m()==0 || A.n()==0) {
    for(int i=0; i<Y.m(); i++)
      Y(i) = beta*Y(i);
  } else {
    iC( dgemv(A.m(), A.n(), alpha, A.data(), X.data(), beta, Y.data()) );
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int dgemv(int m, int n, float alpha, float* A, float* X, float beta, float* Y)
{
  char trans = 'N';
  assert(m!=0 && n!=0);
  int incx = 1;
  int incy = 1;
  sgemv_(&trans, &m, &n, &alpha, A, &m, X, &incx, &beta, Y, &incy);
  return 0;
}

// ---------------------------------------------------------------------- 
int zgemm(cpx alpha, const CpxNumMat& A, const CpxNumMat& B, cpx beta, CpxNumMat& C)
{
  assert( A.m() == C.m() );  assert( A.n() == B.m() );  assert( B.n() == C.n() );
  if(A.m()==0 || A.n()==0 || B.n()==0) { //simply scale
    for(int i=0; i<C.m(); i++)
      for(int j=0; j<C.n(); j++)
	C(i,j) = beta*C(i,j);
  } else {
    iC( zgemm(C.m(), C.n(), A.n(), alpha, A.data(), B.data(), beta, C.data()) );
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int zgemm(int m, int n, int k, cpx alpha, cpx* A, cpx* B, cpx beta, cpx* C)
{
  char transa = 'N';
  char transb = 'N';
  assert(m!=0 && n!=0 && k!=0);
  cgemm_(&transa, &transb, &m, &n, &k,
	 &alpha, A, &m, B, &k, &beta, C, &m);
  return 0;
}
//Y <- a M X + b Y
// ---------------------------------------------------------------------- 
int zgemv(cpx alpha, const CpxNumMat& A, const CpxNumVec& X, cpx beta, CpxNumVec& Y)
{
  assert(Y.m() == A.m());
  assert(A.n() == X.m());
  if(A.m()==0 || A.n()==0) {
    for(int i=0; i<Y.m(); i++)
      Y(i) = beta*Y(i);
  } else {
    iC( zgemv(A.m(), A.n(), alpha, A.data(), X.data(), beta, Y.data()) );
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int zgemv(int m, int n, cpx alpha, cpx* A, cpx* X, cpx beta, cpx* Y)
{
  char trans = 'N';
  int incx = 1;
  int incy = 1;
  assert(m!=0 && n!=0);
  //cerr<<sizeof(int)<<" "<<sizeof(long int)<<endl;
  cgemv_(&trans, &m, &n, &alpha, A, &m, X, &incx, &beta, Y, &incy);
  return 0;
}

// ---------------------------------------------------------------------- 
int dgmres(int (*A)(const FltNumVec&, FltNumVec&), const FltNumVec& b, const FltNumVec& x0,
	   int restart, float tol, int maxit, int print,
	   FltNumVec& x, int& flag, float& relres, int& niter, vector<float>& resvec)
{
  int n = b.m();
  int m = restart;
  
  FltNumMat V(n,m+1);  setvalue(V,0.0f);
  FltNumMat H(m+1,m);  setvalue(H,0.0f);
  
  float bnrm2 = 0;  for(int a=0; a<n; a++)    bnrm2 = bnrm2 + b(a)*b(a);
  bnrm2 = sqrt(bnrm2);
  resvec.clear();
  
  FltNumVec tmp(n);  for(int a=0; a<n; a++)    tmp(a) = 0;
  x = x0;
  float xnrm2 = 0;  for(int a=0; a<n; a++)    xnrm2 = xnrm2 + x(a)*x(a);
  xnrm2 = sqrt(xnrm2);
  if(xnrm2 > 1e-16) {
    iC( (*A)(x,tmp) );
  }
  FltNumVec r(n);  for(int a=0; a<n; a++)    r(a) = b(a) - tmp(a);
  float beta=0;  for(int a=0; a<n; a++)    beta = beta + r(a)*r(a);
  beta = sqrt(beta);
  float res = beta; if(print==1) cerr<<"Iter "<<resvec.size()<<": "<<res<<endl;
  resvec.push_back(res);
  float err = res/bnrm2;
  
  int iter = 0;
  while(1) {
    for(int a=0; a<n; a++)      V(a,0) = r(a)/beta;
    int j = 0;
    FltNumVec y;
    FltNumMat Hj;
    while(1) {
      FltNumVec Vj(n);     for(int a=0; a<n; a++)	Vj(a) = V(a,j);
      FltNumVec w(n);      setvalue(w,0.0f);
      iC( (*A)(Vj,w) );
      for(int k=0; k<=j; k++) {
	float sum = 0;	for(int a=0; a<n; a++)	  sum = sum + V(a,k)*w(a);
	H(k,j) = sum;
	for(int a=0; a<n; a++)	  w(a) = w(a) - sum*V(a,k);
      }
      float nw=0;      for(int a=0; a<n; a++)	nw = nw + w(a)*w(a);
      nw = sqrt(nw);
      H(j+1,j) = nw;
      for(int a=0; a<n; a++)	V(a,j+1) = w(a) / nw;
      FltNumVec be(j+2);      for(int a=0; a<j+2; a++)	be(a) = 0;
      be(0) = beta;
      y.resize(j+1);
      Hj.resize(j+2,j+1);
      for(int a=0; a<j+2; a++)	for(int c=0; c<j+1; c++)	  Hj(a,c) = H(a,c);
      //SOLVE
      {
	int m = j+2;
	int n = j+1;
	int nrhs = 1;
	FltNumMat Hjtmp(Hj);	//cpx* aptr = Hjtmp.data();
	int lda = j+2;
	FltNumVec betmp(be);	//cpx* bptr = betmp.data();
	int ldb = j+2;
	FltNumVec s(j+2);
	float rcond = 0;
	int rank;
	FltNumVec work(10*(j+2));
	int lwork = 10*(j+2);
	//FltNumVec rwork(10*(j+2));
	int info;
	sgelss_(&m,&n,&nrhs,Hjtmp.data(),&lda,betmp.data(),&ldb,s.data(),&rcond,&rank,work.data(),
		&lwork,&info);
	//cerr<<info<<endl;
	for(int a=0; a<j+1; a++)	  y(a) = betmp(a);
      }
      iC( dgemv(-1.0f, Hj, y, 1, be) );
      float res=0;      for(int a=0; a<j+2; a++)	res = res + be(a)*be(a);
      res = sqrt(res);      if(print==1) cerr<<"Iter "<<resvec.size()<<": "<<res<<endl;
      resvec.push_back(res);
      err = res/bnrm2;
      if(err<tol || j==m-1)
	break;
      j=j+1;
    }
    FltNumMat Vj(n,j+1,false,V.data());
    iC( dgemv(1.0f, Vj, y, 1.0, x) );
    FltNumVec tmp(j+2);
    iC( dgemv(1.0, Hj, y, 0.0, tmp) );
    FltNumMat Vj1(n,j+2,false,V.data());
    iC( dgemv(-1.0, Vj1, tmp, 1.0, r) );
    beta = 0;    for(int a=0; a<n; a++)      beta = beta + r(a)*r(a);
    beta = sqrt(beta);
    if(err<tol || iter==maxit-1)
      break;
    iter++;
  }
  flag = (err>tol);
  relres = beta;
  niter = iter+1;
  
  return 0;
}

// ---------------------------------------------------------------------- 
int zgmres(int (*A)(const CpxNumVec&, CpxNumVec&), const CpxNumVec& b, const CpxNumVec& x0,
	   int restart, float tol, int maxit, int print,
	   CpxNumVec& x, int& flag, float& relres, int& niter, vector<float>& resvec)
{
  int n = b.m();
  int m = restart;
  
  CpxNumMat V(n,m+1);  setvalue(V,cpx(0,0));
  CpxNumMat H(m+1,m);  setvalue(H,cpx(0,0));
  
  float bnrm2 = 0;  for(int a=0; a<n; a++)    bnrm2 = bnrm2 + abs(b(a)*b(a));
  bnrm2 = sqrt(bnrm2);
  resvec.clear();
  
  CpxNumVec tmp(n);  for(int a=0; a<n; a++)    tmp(a) = 0;
  x = x0;
  float xnrm2 = 0;  for(int a=0; a<n; a++)    xnrm2 = xnrm2 + abs(x(a)*x(a));
  xnrm2 = sqrt(xnrm2);
  if(xnrm2 > 1e-16) {
    iC( (*A)(x,tmp) );
  }
  CpxNumVec r(n);  for(int a=0; a<n; a++)    r(a) = b(a) - tmp(a);
  float beta=0;  for(int a=0; a<n; a++)    beta = beta + abs(r(a)*r(a));
  beta = sqrt(beta);
  float res = beta; if(print==1) cerr<<"Iter "<<resvec.size()<<": "<<res<<endl;
  resvec.push_back(res);
  float err = res/bnrm2;
  
  int iter = 0;
  while(1) {
    for(int a=0; a<n; a++)      V(a,0) = r(a)/beta;
    int j = 0;
    CpxNumVec y;
    CpxNumMat Hj;
    while(1) {
      CpxNumVec Vj(n);      for(int a=0; a<n; a++)	Vj(a) = V(a,j);
      CpxNumVec w(n);      setvalue(w,cpx(0,0));
      iC( (*A)(Vj,w) );
      for(int k=0; k<=j; k++) {
	cpx sum = 0;	for(int a=0; a<n; a++)	  sum = sum + conj(V(a,k))*w(a);
	H(k,j) = sum;
	for(int a=0; a<n; a++)	  w(a) = w(a) - sum*V(a,k);
      }
      float nw=0;      for(int a=0; a<n; a++)	nw = nw + abs(w(a)*w(a));
      nw = sqrt(nw);
      H(j+1,j) = nw;
      for(int a=0; a<n; a++)	V(a,j+1) = w(a) / nw;
      CpxNumVec be(j+2);      for(int a=0; a<j+2; a++)	be(a) = 0;
      be(0) = beta;
      y.resize(j+1);
      Hj.resize(j+2,j+1);
      for(int a=0; a<j+2; a++)	for(int c=0; c<j+1; c++)	  Hj(a,c) = H(a,c);
      //SOLVE
      {
	int m = j+2;
	int n = j+1;
	int nrhs = 1;
	CpxNumMat Hjtmp(Hj);	//cpx* aptr = Hjtmp.data();
	int lda = j+2;
	CpxNumVec betmp(be);	//cpx* bptr = betmp.data();
	int ldb = j+2;
	FltNumVec s(j+2);
	float rcond = 0;
	int rank;
	CpxNumVec work(10*(j+2));
	int lwork = 10*(j+2);
	FltNumVec rwork(10*(j+2));
	int info;
	cgelss_(&m,&n,&nrhs,
		(MKL_Complex8*) Hjtmp.data(),&lda,
		(MKL_Complex8*) betmp.data(),&ldb,
		s.data(),&rcond,&rank,
		(MKL_Complex8*) work.data(),&lwork,
		rwork.data(),&info);
	for(int a=0; a<j+1; a++)	  y(a) = betmp(a);
      }
      iC( zgemv(-1.0, Hj, y, 1, be) );
      float res=0;      for(int a=0; a<j+2; a++)	res = res + abs(be(a)*be(a));
      res = sqrt(res);      if(print==1) cerr<<"Iter "<<resvec.size()<<": "<<res<<endl;
      resvec.push_back(res);
      err = res/bnrm2;
      if(err<tol || j==m-1)
	break;
      j=j+1;
    }
    CpxNumMat Vj(n,j+1,false,V.data());
    iC( zgemv(1.0, Vj, y, 1.0, x) );
    CpxNumVec tmp(j+2);
    iC( zgemv(1.0, Hj, y, 0.0, tmp) );
    CpxNumMat Vj1(n,j+2,false,V.data());
    iC( zgemv(-1.0, Vj1, tmp, 1.0, r) );
    beta = 0;    for(int a=0; a<n; a++)      beta = beta + abs(r(a)*r(a));
    beta = sqrt(beta);
    if(err<tol || iter==maxit-1)
      break;
    iter++;
  }
  flag = (err>tol);
  relres = beta;
  niter = iter+1;
  
  return 0;
}

// ---------------------------------------------------------------------- 
int pinv(const FltNumMat& M, float eps, FltNumMat& R)
{
  //svd
  int m = M.m();
  int n = M.n();
  int k = min(m,n);
  
  if(m==0 || n==0) { //LEXING: IMPORTANT
    R.resize(n,m);    setvalue(R, 0.0f);
  } else {
    FltNumMat U(m,k);
    FltNumVec S(k);
    FltNumMat VT(k,n);
    {
      char jobu  = 'S';
      char jobvt = 'S';
      FltNumMat MC(M);
      int lwork = 20*max(m,n);
      FltNumVec work(lwork);
      int info;
      sgesvd_(&jobu, &jobvt, &m, &n, MC.data(), &m, S.data(), U.data(), &m, VT.data(), &k, work.data(), &lwork, &info);    iA(info==0);
    }
    //threshold
    float cutoff=eps*S(0); //relative thresholding
    int p = 0;
    for(int a=0; a<k; a++)    if(abs(S(a))>cutoff)      p++;
    
    //switch
    FltNumMat UT(p,m);
    for(int i=0; i<m; i++)
      for(int j=0; j<p; j++)
	UT(j,i) = U(i,j);
    FltNumMat V(n,p);
    for(int i=0; i<p; i++)
      for(int j=0; j<n; j++)
	V(j,i) = VT(i,j);
    FltNumVec IS(p);
    for(int i=0; i<p; i++)
      IS(i) = 1.0/S(i);
    
    //multiply back
    for(int i=0; i<n; i++)
      for(int j=0; j<p; j++)
	V(i,j) *= IS(j);
    R.resize(n,m);
    iC( dgemm(1.0, V, UT, 0.0, R) );
  }
  return 0;
}

// complex version
int pinv(const CpxNumMat& M, float eps, CpxNumMat& R)
{
  //svd
  int m = M.m();
  int n = M.n();
  int k = min(m,n);
  
  if(m==0 || n==0) { //LEXING: IMPORTANT
    R.resize(n,m);    setvalue(R, cpx(0,0));
  } else {
    CpxNumMat U(m,k);
    FltNumVec S(k);
    CpxNumMat VT(k,n);
    {
      char jobu  = 'S';
      char jobvt = 'S';
      CpxNumMat MC(M);
      int lwork = 20*max(m,n);
      CpxNumVec work(lwork);
      FltNumVec rwork(lwork);
      int info;
      cgesvd_(&jobu, &jobvt, &m, &n, 
	      (MKL_Complex8*) MC.data(), &m, S.data(), 
	      (MKL_Complex8*) U.data(), &m, 
	      (MKL_Complex8*) VT.data(), &k, 
	      (MKL_Complex8*) work.data(), &lwork, 
	      rwork.data(), &info);    
      iA(info==0);
    }
    //threshold
    float cutoff=eps*S(0); //relative thresholding
    int p = 0;
    for(int a=0; a<k; a++)    if(abs(S(a))>cutoff)      p++;
    
    //switch
    CpxNumMat UT(p,m);
    for(int i=0; i<m; i++)
      for(int j=0; j<p; j++)
	UT(j,i) = conj(U(i,j));
    CpxNumMat V(n,p);
    for(int i=0; i<p; i++)
      for(int j=0; j<n; j++)
	V(j,i) = conj(VT(i,j));
    FltNumVec IS(p);
    for(int i=0; i<p; i++)
      IS(i) = 1.0/S(i);
    
    //multiply back
    for(int i=0; i<n; i++)
      for(int j=0; j<p; j++)
	V(i,j) *= IS(j);
    R.resize(n,m);
    iC( zgemm(1.0, V, UT, 0.0, R) );
  }
  return 0;
}

// ---------------------------------------------------------------------- 
int lowrank(int m, int n, int (*sample)(vector<int>&, vector<int>&, FltNumMat&), float eps, int npk, 
	    vector<int>& cidx, vector<int>& ridx, FltNumMat& mid)
{
  iA(m>0 && n>0);
  {
    int nc = min(npk,n);
    vector<int> cs(nc);
    for(int k=0; k<nc; k++)      cs[k] = int( floor(drand48()*(n-nc)) );
    sort(cs.begin(), cs.end());
    for(int k=0; k<nc; k++)      cs[k] += k;
    //
    //for(int k=0; k<cidx.size(); k++)      cs.push_back(cidx[k]);
    //sort(cs.begin(), cs.end());
    //vector<int>::iterator newend = unique(cs.begin(), cs.end());
    //cs.resize(newend-cs.begin());
    //
    vector<int> rs(m);
    for(int k=0; k<m; k++)      rs[k] = k;
    //
    FltNumMat M2tmp;    iC( (*sample)(rs, cs, M2tmp) );
    FltNumMat M2(M2tmp.n(), M2tmp.m());
    for(int i=0; i<M2tmp.m(); i++)      for(int j=0; j<M2tmp.n(); j++)	M2(j,i) = M2tmp(i,j);
    //
    int m = M2.m();
    int n = M2.n();
    int lda = m;
    NumVec<int> jpvt(n);    setvalue(jpvt, int(0));
    FltNumVec tau(max(m,n));
    FltNumVec work(3*n);
    int lwork = 3*n;
    int info;
    sgeqp3_(&m, &n, M2.data(), &lda, jpvt.data(), tau.data(), work.data(), &lwork, &info);    iA(info==0);
    float cutoff = eps*abs(M2(0,0));
    int cnt=0;
    for(int k=0; k<min(m,n); k++)      if(abs(M2(k,k))>cutoff)	cnt++;
    ridx.resize(cnt);
    for(int k=0; k<cnt; k++)      ridx[k] = (jpvt(k)-1);
    cerr<<"ROWS "; for(int k=0; k<cnt; k++)      cerr<<ridx[k]<<" ";    cerr<<endl;
  }
  {
    int nr = min(npk,m);
    vector<int> rs(nr);
    for(int k=0; k<nr; k++)      rs[k] = int( floor(drand48()*(m-nr)) );
    sort(rs.begin(), rs.end());
    for(int k=0; k<nr; k++)      rs[k] += k;
    //
    for(unsigned int k=0; k<ridx.size(); k++)      rs.push_back(ridx[k]);
    sort(rs.begin(), rs.end());
    vector<int>::iterator newend = unique(rs.begin(), rs.end());
    rs.resize(newend-rs.begin());
    //
    vector<int> cs(n);
    for(int k=0; k<n; k++)      cs[k] = k;
    //
    FltNumMat M1;    iC( (*sample)(rs, cs, M1) );
    //
    int m = M1.m();
    int n = M1.n();
    int lda = m;
    NumVec<int> jpvt(n);      setvalue(jpvt, int(0));
    FltNumVec tau(max(m,n));
    FltNumVec work(3*n);
    int lwork = 3*n;
    int info;
    sgeqp3_(&m, &n, M1.data(), &lda, jpvt.data(), tau.data(), work.data(), &lwork, &info);    iA(info==0);
    float cutoff = eps*abs(M1(0,0)); //the diagonal element
    int cnt=0;
    for(int k=0; k<min(m,n); k++)	if(abs(M1(k,k))>cutoff)	  cnt++;
    cidx.resize(cnt);
    for(int k=0; k<cnt; k++)	cidx[k] = (jpvt(k)-1);
    cerr<<"COLS "; for(int k=0; k<cnt; k++)      cerr<<cidx[k]<<" ";    cerr<<endl;
  }
  {
    int nc = min(npk,n);
    vector<int> cs(nc);
    for(int k=0; k<nc; k++)      cs[k] = int( floor(drand48()*(n-nc)) );
    sort(cs.begin(), cs.end());
    for(int k=0; k<nc; k++)      cs[k] += k;
    for(unsigned int k=0; k<cidx.size(); k++)      cs.push_back(cidx[k]);
    sort(cs.begin(), cs.end());
    vector<int>::iterator csnewend = unique(cs.begin(), cs.end());
    cs.resize(csnewend-cs.begin());
    //
    int nr = min(npk,m);
    vector<int> rs(nr);
    for(int k=0; k<nr; k++)      rs[k] = int( floor(drand48()*(m-nr)) );
    sort(rs.begin(), rs.end());
    for(int k=0; k<nr; k++)      rs[k] += k;
    for(unsigned int k=0; k<ridx.size(); k++)      rs.push_back(ridx[k]);
    sort(rs.begin(), rs.end());
    vector<int>::iterator rsnewend = unique(rs.begin(), rs.end());
    rs.resize(rsnewend-rs.begin());
    //
    FltNumMat M1;    iC( (*sample)(rs,cidx,M1) );
    FltNumMat IM1;    iC( pinv(M1, (float) 1e-7, IM1) );
    FltNumMat M2;    iC( (*sample)(ridx,cs,M2) );
    FltNumMat IM2;    iC( pinv(M2, (float) 1e-7, IM2) );
    FltNumMat M3;    iC( (*sample)(rs,cs,M3) );
    FltNumMat tmp(M3.m(), IM2.n());
    iC( dgemm(1.0, M3, IM2, 0.0, tmp) );
    mid.resize(IM1.m(), tmp.n());
    iC( dgemm(1.0, IM1, tmp, 0.0, mid) );
  }
  if(1) {
    int nc = min(npk,n);
    vector<int> cs(nc);
    for(int k=0; k<nc; k++)      cs[k] = int( floor(drand48()*n) );
    int nr = min(npk,m);
    vector<int> rs(nr);
    for(int k=0; k<nr; k++)      rs[k] = int( floor(drand48()*m) );
    FltNumMat M1;
    iC( (*sample)(rs,cidx,M1) );
    FltNumMat M2;
    iC( (*sample)(ridx,cs,M2) );
    FltNumMat Mext;
    iC( (*sample)(rs,cs,Mext) );
    FltNumMat Mapp(rs.size(), cs.size());
    FltNumMat tmp(mid.m(),M2.n());
    iC( dgemm(1.0, mid, M2, 0.0, tmp) );
    iC( dgemm(1.0, M1, tmp, 0.0, Mapp) );
    FltNumMat Merr(Mext.m(), Mext.n());
    for(int a=0; a<Mext.m(); a++)
      for(int b=0; b<Mext.n(); b++)
	Merr(a,b) = Mext(a,b) - Mapp(a,b);
    cerr<<"rel err "<<sqrt(energy(Merr))/sqrt(energy(Mext))<<endl;
  }
  return 0;
}

//complex version
int lowrank(int m, int n, int (*sample)(vector<int>&, vector<int>&, CpxNumMat&), float eps, int npk, 
	    vector<int>& cidx, vector<int>& ridx, CpxNumMat& mid)
{
  iA(m>0 && n>0);
  {
    int nc = min(npk,n);
    vector<int> cs(nc);
    for(int k=0; k<nc; k++)      cs[k] = int( floor(drand48()*(n-nc)) );
    sort(cs.begin(), cs.end());
    for(int k=0; k<nc; k++)      cs[k] += k;
    //
    //for(int k=0; k<cidx.size(); k++)      cs.push_back(cidx[k]);
    //sort(cs.begin(), cs.end());
    //vector<int>::iterator newend = unique(cs.begin(), cs.end());
    //cs.resize(newend-cs.begin());
    //
    vector<int> rs(m);
    for(int k=0; k<m; k++)      rs[k] = k;
    //
    CpxNumMat M2tmp;    iC( (*sample)(rs, cs, M2tmp) );
    CpxNumMat M2(M2tmp.n(), M2tmp.m());
    for(int i=0; i<M2tmp.m(); i++)      for(int j=0; j<M2tmp.n(); j++)	M2(j,i) = conj(M2tmp(i,j));
    //
    int m = M2.m();
    int n = M2.n();
    int lda = m;
    NumVec<int> jpvt(n);    setvalue(jpvt, int(0));
    CpxNumVec tau(max(m,n));
    CpxNumVec work(3*n);
    FltNumVec rwork(6*n);
    int lwork = 3*n;
    int info;
    cgeqp3_(&m, &n, 
	    (MKL_Complex8*) M2.data(), &lda, jpvt.data(), 
	    (MKL_Complex8*) tau.data(), 
	    (MKL_Complex8*) work.data(), &lwork, rwork.data(), &info);    
    iA(info==0);
    float cutoff = eps*abs(M2(0,0));
    int cnt=0;
    for(int k=0; k<min(m,n); k++)      if(abs(M2(k,k))>cutoff)	cnt++;
    ridx.resize(cnt);
    for(int k=0; k<cnt; k++)      ridx[k] = (jpvt(k)-1);
    cerr<<"ROWS "; for(int k=0; k<cnt; k++)      cerr<<ridx[k]<<" ";    cerr<<endl;
  }
  {
    int nr = min(npk,m);
    vector<int> rs(nr);
    for(int k=0; k<nr; k++)      rs[k] = int( floor(drand48()*(m-nr)) );
    sort(rs.begin(), rs.end());
    for(int k=0; k<nr; k++)      rs[k] += k;
    //
    for(unsigned int k=0; k<ridx.size(); k++)      rs.push_back(ridx[k]);
    sort(rs.begin(), rs.end());
    vector<int>::iterator newend = unique(rs.begin(), rs.end());
    rs.resize(newend-rs.begin());
    //
    vector<int> cs(n);
    for(int k=0; k<n; k++)      cs[k] = k;
    //
    CpxNumMat M1;    iC( (*sample)(rs, cs, M1) );
    //
    int m = M1.m();
    int n = M1.n();
    int lda = m;
    NumVec<int> jpvt(n);      setvalue(jpvt, int(0));
    CpxNumVec tau(max(m,n));
    CpxNumVec work(3*n);
    FltNumVec rwork(6*n);
    int lwork = 3*n;
    int info;
    cgeqp3_(&m, &n, 
	    (MKL_Complex8*) M1.data(), &lda, jpvt.data(), 
	    (MKL_Complex8*) tau.data(), 
	    (MKL_Complex8*) work.data(), &lwork, rwork.data(), &info);    
    iA(info==0);
    float cutoff = eps*abs(M1(0,0)); //the diagonal element
    int cnt=0;
    for(int k=0; k<min(m,n); k++)	if(abs(M1(k,k))>cutoff)	  cnt++;
    cidx.resize(cnt);
    for(int k=0; k<cnt; k++)	cidx[k] = (jpvt(k)-1);
    cerr<<"COLS "; for(int k=0; k<cnt; k++)      cerr<<cidx[k]<<" ";    cerr<<endl;
  }
  {
    int nc = min(npk,n);
    vector<int> cs(nc);
    for(int k=0; k<nc; k++)      cs[k] = int( floor(drand48()*(n-nc)) );
    sort(cs.begin(), cs.end());
    for(int k=0; k<nc; k++)      cs[k] += k;
    for(unsigned int k=0; k<cidx.size(); k++)      cs.push_back(cidx[k]);
    sort(cs.begin(), cs.end());
    vector<int>::iterator csnewend = unique(cs.begin(), cs.end());
    cs.resize(csnewend-cs.begin());
    //
    int nr = min(npk,m);
    vector<int> rs(nr);
    for(int k=0; k<nr; k++)      rs[k] = int( floor(drand48()*(m-nr)) );
    sort(rs.begin(), rs.end());
    for(int k=0; k<nr; k++)      rs[k] += k;
    for(unsigned int k=0; k<ridx.size(); k++)      rs.push_back(ridx[k]);
    sort(rs.begin(), rs.end());
    vector<int>::iterator rsnewend = unique(rs.begin(), rs.end());
    rs.resize(rsnewend-rs.begin());
    //
    CpxNumMat M1;    iC( (*sample)(rs,cidx,M1) );
    CpxNumMat IM1;    iC( pinv(M1, (float) 1e-7, IM1) );
    CpxNumMat M2;    iC( (*sample)(ridx,cs,M2) );
    CpxNumMat IM2;    iC( pinv(M2, (float) 1e-7, IM2) );
    CpxNumMat M3;    iC( (*sample)(rs,cs,M3) );
    CpxNumMat tmp(M3.m(), IM2.n());
    iC( zgemm(1.0, M3, IM2, 0.0, tmp) );
    mid.resize(IM1.m(), tmp.n());
    iC( zgemm(1.0, IM1, tmp, 0.0, mid) );
  }
  if(1) {
    int nc = min(npk,n);
    vector<int> cs(nc);
    for(int k=0; k<nc; k++)      cs[k] = int( floor(drand48()*n) );
    int nr = min(npk,m);
    vector<int> rs(nr);
    for(int k=0; k<nr; k++)      rs[k] = int( floor(drand48()*m) );
    CpxNumMat M1;
    iC( (*sample)(rs,cidx,M1) );
    CpxNumMat M2;
    iC( (*sample)(ridx,cs,M2) );
    CpxNumMat Mext;
    iC( (*sample)(rs,cs,Mext) );
    CpxNumMat Mapp(rs.size(), cs.size());
    CpxNumMat tmp(mid.m(),M2.n());
    iC( zgemm(1.0, mid, M2, 0.0, tmp) );
    iC( zgemm(1.0, M1, tmp, 0.0, Mapp) );
    CpxNumMat Merr(Mext.m(), Mext.n());
    for(int a=0; a<Mext.m(); a++)
      for(int b=0; b<Mext.n(); b++)
	Merr(a,b) = Mext(a,b) - Mapp(a,b);
    cerr<<"lowrank rel err "<<sqrt(energy(Merr))/sqrt(energy(Mext))<<endl;
  }
  return 0;
}


//-------------------------------------------------------------------------
//just the transpose,
int ztran(const CpxNumMat& A, CpxNumMat& B)
{
  B.resize(A.n(), A.m());
  for(int i=0; i<B.m(); i++)
    for(int j=0; j<B.n(); j++)
      B(i,j) = A(j,i);
  return 0;
}

//-------------------------------------------------------------------------
//shift left
int shiftleft(const CpxNumTns& A, CpxNumTns& B)
{
  iA( B.m()==A.n() && B.n()==A.p() && B.p()==A.m() );  //B.resize(A.n(), A.p(), A.m());
  for(int i=0; i<B.m(); i++)
    for(int j=0; j<B.n(); j++)
      for(int k=0; k<B.p(); k++)
	B(i,j,k) = A(k,i,j);
  return 0;
}


//---------------------------------------------------------------------------------------------------
//////////Double



// ---------------------------------------------------------------------- 
int ddgemm(double alpha, const DblNumMat& A, const DblNumMat& B, double beta, DblNumMat& C)
{
  assert( A.m() == C.m() );  assert( A.n() == B.m() );  assert( B.n() == C.n() );
  if(A.m()==0 || A.n()==0 || B.n()==0) { //simply scale
    for(int i=0; i<C.m(); i++)
      for(int j=0; j<C.n(); j++)
        C(i,j) = beta*C(i,j);
  } else {
    iC( ddgemm(C.m(), C.n(), A.n(), alpha, A.data(), B.data(), beta, C.data()) );
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int ddgemm(int m, int n, int k, double alpha, double* A, double* B, double beta, double* C)
{
  char transa = 'N';
  char transb = 'N';
  assert(m!=0 && n!=0 && k!=0);
  dgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &m, B, &k, &beta, C, &m);
  return 0;
}
//Y <- a M X + b Y
// ---------------------------------------------------------------------- 

int ddgemv(double alpha, const DblNumMat& A, const DblNumVec& X, double beta, DblNumVec& Y)
{
  assert(Y.m() == A.m());
  assert(A.n() == X.m());
  if(A.m()==0 || A.n()==0) {
    for(int i=0; i<Y.m(); i++)
      Y(i) = beta*Y(i);
  } else {
    iC( ddgemv(A.m(), A.n(), alpha, A.data(), X.data(), beta, Y.data()) );
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int ddgemv(int m, int n, double alpha, double* A, double* X, double beta, double* Y)
{
  char trans = 'N';
  assert(m!=0 && n!=0);
  int incx = 1;
  int incy = 1;
  dgemv_(&trans, &m, &n, &alpha, A, &m, X, &incx, &beta, Y, &incy);
  return 0;
}

// double complex version
// ---------------------------------------------------------------------- 
int zzgemm(zpx alpha, const ZpxNumMat& A, const ZpxNumMat& B, zpx beta, ZpxNumMat& C)
{
  assert( A.m() == C.m() );  assert( A.n() == B.m() );  assert( B.n() == C.n() );
  if(A.m()==0 || A.n()==0 || B.n()==0) { //simply scale
    for(int i=0; i<C.m(); i++)
      for(int j=0; j<C.n(); j++)
	C(i,j) = beta*C(i,j);
  } else {
    iC( zzgemm(C.m(), C.n(), A.n(), alpha, A.data(), B.data(), beta, C.data()) );
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int zzgemm(int m, int n, int k, zpx alpha, zpx* A, zpx* B, zpx beta, zpx* C)
{
  char transa = 'N';
  char transb = 'N';
  assert(m!=0 && n!=0 && k!=0);
  zgemm_(&transa, &transb, &m, &n, &k,
	 &alpha, A, &m, B, &k, &beta, C, &m);
  return 0;
}
//Y <- a M X + b Y
// ---------------------------------------------------------------------- 

int zzgemv(zpx alpha, const ZpxNumMat& A, const ZpxNumVec& X, zpx beta, ZpxNumVec& Y)
{
  assert(Y.m() == A.m());
  assert(A.n() == X.m());
  if(A.m()==0 || A.n()==0) {
    for(int i=0; i<Y.m(); i++)
      Y(i) = beta*Y(i);
  } else {
    iC( zzgemv(A.m(), A.n(), alpha, A.data(), X.data(), beta, Y.data()) );
  }
  return 0;
}
// ---------------------------------------------------------------------- 
int zzgemv(int m, int n, zpx alpha, zpx* A, zpx* X, zpx beta, zpx* Y)
{
  char trans = 'N';
  int incx = 1;
  int incy = 1;
  assert(m!=0 && n!=0);
  //cerr<<sizeof(int)<<" "<<sizeof(long int)<<endl;
  zgemv_(&trans, &m, &n, &alpha, A, &m, X, &incx, &beta, Y, &incy);
  return 0;
}

// ---------------------------------------------------------------------- 
int ddgmres(int (*A)(const DblNumVec&, DblNumVec&), const DblNumVec& b, const DblNumVec& x0,
           int restart, double tol, int maxit, int print,
           DblNumVec& x, int& flag, double& relres, int& niter, vector<double>& resvec)
{
  int n = b.m();
  int m = restart;

  DblNumMat V(n,m+1);  setvalue(V,0.0);
  DblNumMat H(m+1,m);  setvalue(H,0.0);

  double bnrm2 = 0;  for(int a=0; a<n; a++)    bnrm2 = bnrm2 + b(a)*b(a);
  bnrm2 = sqrt(bnrm2);
  resvec.clear();

  DblNumVec tmp(n);  for(int a=0; a<n; a++)    tmp(a) = 0;
  x = x0;
  double xnrm2 = 0;  for(int a=0; a<n; a++)    xnrm2 = xnrm2 + x(a)*x(a);
  xnrm2 = sqrt(xnrm2);
  if(xnrm2 > 1e-16) {
    iC( (*A)(x,tmp) );
  }
  DblNumVec r(n);  for(int a=0; a<n; a++)    r(a) = b(a) - tmp(a);
  double beta=0;  for(int a=0; a<n; a++)    beta = beta + r(a)*r(a);
  beta = sqrt(beta);
  double res = beta; if(print==1) cerr<<"Iter "<<resvec.size()<<": "<<res<<endl;
  resvec.push_back(res);
  double err = res/bnrm2;

  int iter = 0;
  while(1) {
    for(int a=0; a<n; a++)      V(a,0) = r(a)/beta;
    int j = 0;
    DblNumVec y;
    DblNumMat Hj;
    while(1) {
      DblNumVec Vj(n);      for(int a=0; a<n; a++)      Vj(a) = V(a,j);
      DblNumVec w(n);      setvalue(w,0.0);
      iC( (*A)(Vj,w) );
      for(int k=0; k<=j; k++) {
        double sum = 0; for(int a=0; a<n; a++)    sum = sum + V(a,k)*w(a);
        H(k,j) = sum;
        for(int a=0; a<n; a++)    w(a) = w(a) - sum*V(a,k);
      }
      double nw=0;      for(int a=0; a<n; a++)  nw = nw + w(a)*w(a);
      nw = sqrt(nw);
      H(j+1,j) = nw;
      for(int a=0; a<n; a++)    V(a,j+1) = w(a) / nw;
      DblNumVec be(j+2);      for(int a=0; a<j+2; a++)  be(a) = 0;
      be(0) = beta;
      y.resize(j+1);
      Hj.resize(j+2,j+1);
      for(int a=0; a<j+2; a++)  for(int c=0; c<j+1; c++)          Hj(a,c) = H(a,c);
      //SOLVE
      {
        int m = j+2;
        int n = j+1;
        int nrhs = 1;
        DblNumMat Hjtmp(Hj);    //cpx* aptr = Hjtmp.data();
        int lda = j+2;
        DblNumVec betmp(be);    //cpx* bptr = betmp.data();
        int ldb = j+2;
        DblNumVec s(j+2);
        double rcond = 0;
        int rank;
        DblNumVec work(10*(j+2));
        int lwork = 10*(j+2);
        //DblNumVec rwork(10*(j+2));
        int info;
        dgelss_(&m,&n,&nrhs,Hjtmp.data(),&lda,betmp.data(),&ldb,s.data(),&rcond,&rank,work.data(),
                &lwork,&info);
        //cerr<<info<<endl;
        for(int a=0; a<j+1; a++)          y(a) = betmp(a);
      }
      iC( ddgemv(-1.0, Hj, y, 1, be) );
      double res=0;      for(int a=0; a<j+2; a++)       res = res + be(a)*be(a);
      res = sqrt(res);      if(print==1) cerr<<"Iter "<<resvec.size()<<": "<<res<<endl;
      resvec.push_back(res);
      err = res/bnrm2;
      if(err<tol || j==m-1)
        break;
      j=j+1;
    }
    DblNumMat Vj(n,j+1,false,V.data());
    iC( ddgemv(1.0, Vj, y, 1.0, x) );
    DblNumVec tmp(j+2);
    iC( ddgemv(1.0, Hj, y, 0.0, tmp) );
    DblNumMat Vj1(n,j+2,false,V.data());
    iC( ddgemv(-1.0, Vj1, tmp, 1.0, r) );
    beta = 0;    for(int a=0; a<n; a++)      beta = beta + r(a)*r(a);
    beta = sqrt(beta);
    if(err<tol || iter==maxit-1)
      break;
    iter++;
  }
  flag = (err>tol);
  relres = beta;
  niter = iter+1;

  return 0;
}
// ---------------------------------------------------------------------- 
int ddpinv(const DblNumMat& M, double eps, DblNumMat& R)
{
  //svd
  int m = M.m();
  int n = M.n();
  int k = min(m,n);

  if(m==0 || n==0) { //LEXING: IMPORTANT
    R.resize(n,m);    setvalue(R, 0.0);
  } else {
    DblNumMat U(m,k);
    DblNumVec S(k);
    DblNumMat VT(k,n);
    {
      char jobu  = 'S';
      char jobvt = 'S';
      DblNumMat MC(M);
      int lwork = 20*max(m,n);
      DblNumVec work(lwork);
      int info;
      dgesvd_(&jobu, &jobvt, &m, &n, MC.data(), &m, S.data(), U.data(), &m, VT.data(), &k, work.data(), &lwork, &info);    iA(info==0);
    }
    //threshold
    double cutoff=eps*S(0); //relative thresholding
    int p = 0;
    for(int a=0; a<k; a++)    if(abs(S(a))>cutoff)      p++;

    //switch
    DblNumMat UT(p,m);
    for(int i=0; i<m; i++)
      for(int j=0; j<p; j++)
        UT(j,i) = U(i,j);

    DblNumMat V(n,p);
    for(int i=0; i<p; i++)
      for(int j=0; j<n; j++)
        V(j,i) = VT(i,j);
    DblNumVec IS(p);
    for(int i=0; i<p; i++)
      IS(i) = 1.0/S(i);

    //multiply back
    for(int i=0; i<n; i++)
      for(int j=0; j<p; j++)
        V(i,j) *= IS(j);
    R.resize(n,m);
    iC( ddgemm(1.0, V, UT, 0.0, R) );
  }
  return 0;
}
//complex version
// ---------------------------------------------------------------------- 
int ddpinv(const ZpxNumMat& M, double eps, ZpxNumMat& R)
{
  //svd
  int m = M.m();
  int n = M.n();
  int k = min(m,n);

  if(m==0 || n==0) { //LEXING: IMPORTANT
    R.resize(n,m);    setvalue(R, zpx(0.,0.)); //Junzhe: important!!! zpx(0.,0.)?
  } else {
    ZpxNumMat U(m,k);
    DblNumVec S(k); //Junzhe: important!!! singular values are real
    ZpxNumMat VT(k,n);
    {
      char jobu  = 'S';
      char jobvt = 'S';
      ZpxNumMat MC(M);
      int lwork = 20*max(m,n);
      ZpxNumVec work(lwork);
      DblNumVec rwork(lwork);
      int info;
      zgesvd_(&jobu, &jobvt, &m, &n, 
	      (MKL_Complex16*) MC.data(), &m, S.data(), 
	      (MKL_Complex16*) U.data(), &m, 
	      (MKL_Complex16*) VT.data(), &k, 
	      (MKL_Complex16*) work.data(), &lwork, 
	      rwork.data(), &info);    
      iA(info==0);
    }
    //threshold
    double cutoff=eps*S(0); //relative thresholding
    int p = 0;
    for(int a=0; a<k; a++)    if(abs(S(a))>cutoff)      p++;

    //switch
    ZpxNumMat UT(p,m);
    for(int i=0; i<m; i++)
      for(int j=0; j<p; j++)
	UT(j,i) = conj(U(i,j));
    ZpxNumMat V(n,p);
    for(int i=0; i<p; i++)
      for(int j=0; j<n; j++)
	V(j,i) = conj(VT(i,j));
    DblNumVec IS(p);
    for(int i=0; i<p; i++)
      IS(i) = 1.0/S(i);
    
    //multiply back
    for(int i=0; i<n; i++)
      for(int j=0; j<p; j++)
	V(i,j) *= IS(j);
    R.resize(n,m);
    iC( zzgemm(1.0, V, UT, 0.0, R) );
  }
  return 0;
}

// ---------------------------------------------------------------------- 
int ddlowrank(int m, int n, int (*sample)(vector<int>&, vector<int>&, DblNumMat&), double eps, int npk,
            vector<int>& cidx, vector<int>& ridx, DblNumMat& mid)
{
  iA(m>0 && n>0);
  {
    int nc = min(npk,n);
    vector<int> cs(nc);
    for(int k=0; k<nc; k++)      cs[k] = int( floor(drand48()*(n-nc)) );
    sort(cs.begin(), cs.end());
    for(int k=0; k<nc; k++)      cs[k] += k;
    vector<int> rs(m);
    for(int k=0; k<m; k++)      rs[k] = k;
    DblNumMat M2tmp;    iC( (*sample)(rs, cs, M2tmp) );
    DblNumMat M2(M2tmp.n(), M2tmp.m());
    for(int i=0; i<M2tmp.m(); i++)      for(int j=0; j<M2tmp.n(); j++)  M2(j,i) = M2tmp(i,j);
    int m = M2.m();
    int n = M2.n();
    int lda = m;
    NumVec<int> jpvt(n);    setvalue(jpvt, int(0));
    DblNumVec tau(max(m,n));
    DblNumVec work(3*n);
    int lwork = 6*n;
    int info;
    dgeqp3_(&m, &n, M2.data(), &lda, jpvt.data(), tau.data(), work.data(), &lwork, &info);    iA(info==0);
    double cutoff = eps*abs(M2(0,0));
    int cnt=0;
    for(int k=0; k<min(m,n); k++)      if(abs(M2(k,k))>cutoff)  cnt++;
    ridx.resize(cnt);
    for(int k=0; k<cnt; k++)      ridx[k] = (jpvt(k)-1);
    cerr<<"ROWS "; for(int k=0; k<cnt; k++)      cerr<<ridx[k]<<" ";    cerr<<endl;
  }
  {
    int nr = min(npk,m);
    vector<int> rs(nr);
    for(int k=0; k<nr; k++)      rs[k] = int( floor(drand48()*(m-nr)) );
    sort(rs.begin(), rs.end());
    for(int k=0; k<nr; k++)      rs[k] += k;
    for(unsigned int k=0; k<ridx.size(); k++)      rs.push_back(ridx[k]);
    sort(rs.begin(), rs.end());
    vector<int>::iterator newend = unique(rs.begin(), rs.end());
    rs.resize(newend-rs.begin());
    vector<int> cs(n);
    for(int k=0; k<n; k++)      cs[k] = k;
    DblNumMat M1;    iC( (*sample)(rs, cs, M1) );
    int m = M1.m();
    int n = M1.n();
    int lda = m;
    NumVec<int> jpvt(n);      setvalue(jpvt, int(0));
    DblNumVec tau(max(m,n));
    DblNumVec work(3*n);
    int lwork = 6*n;
    int info;
    dgeqp3_(&m, &n, M1.data(), &lda, jpvt.data(), tau.data(), work.data(), &lwork, &info);    iA(info==0);
    double cutoff = eps*abs(M1(0,0));    int cnt=0;
    for(int k=0; k<min(m,n); k++)       if(abs(M1(k,k))>cutoff)   cnt++;
    cidx.resize(cnt);
    for(int k=0; k<cnt; k++)    cidx[k] = (jpvt(k)-1);
    cerr<<"COLS "; for(int k=0; k<cnt; k++)      cerr<<cidx[k]<<" ";    cerr<<endl;
  }

  {
    int nc = min(npk,n);
    vector<int> cs(nc);
    for(int k=0; k<nc; k++)      cs[k] = int( floor(drand48()*(n-nc)) );
    sort(cs.begin(), cs.end());
    for(int k=0; k<nc; k++)      cs[k] += k;
    for(unsigned int k=0; k<cidx.size(); k++)      cs.push_back(cidx[k]);
    sort(cs.begin(), cs.end());
    vector<int>::iterator csnewend = unique(cs.begin(), cs.end());
    cs.resize(csnewend-cs.begin());
    int nr = min(npk,m);
    vector<int> rs(nr);
    for(int k=0; k<nr; k++)      rs[k] = int( floor(drand48()*(m-nr)) );
    sort(rs.begin(), rs.end());
    for(int k=0; k<nr; k++)      rs[k] += k;
    for(unsigned int k=0; k<ridx.size(); k++)      rs.push_back(ridx[k]);
    sort(rs.begin(), rs.end());
    vector<int>::iterator rsnewend = unique(rs.begin(), rs.end());
    rs.resize(rsnewend-rs.begin());
    DblNumMat M1;    iC( (*sample)(rs,cidx,M1) );
    DblNumMat IM1;    iC( ddpinv(M1, 1e-16, IM1) );
    DblNumMat M2;    iC( (*sample)(ridx,cs,M2) );
    DblNumMat IM2;    iC( ddpinv(M2, 1e-16, IM2) );
    DblNumMat M3;    iC( (*sample)(rs,cs,M3) );
    DblNumMat tmp(M3.m(), IM2.n());
    iC( ddgemm(1.0, M3, IM2, 0.0, tmp) );
    mid.resize(IM1.m(), tmp.n());
    iC( ddgemm(1.0, IM1, tmp, 0.0, mid) );
  }
  if(1) {
    int nc = min(npk,n);
    vector<int> cs(nc);
    for(int k=0; k<nc; k++)      cs[k] = int( floor(drand48()*n) );
    int nr = min(npk,m);
    vector<int> rs(nr);
    for(int k=0; k<nr; k++)      rs[k] = int( floor(drand48()*m) );
    DblNumMat M1;
    iC( (*sample)(rs,cidx,M1) );
    DblNumMat M2;
    iC( (*sample)(ridx,cs,M2) );
    DblNumMat Mext;
    iC( (*sample)(rs,cs,Mext) );
    DblNumMat Mapp(rs.size(), cs.size());
    DblNumMat tmp(mid.m(),M2.n());
    iC( ddgemm(1.0, mid, M2, 0.0, tmp) );
    iC( ddgemm(1.0, M1, tmp, 0.0, Mapp) );
    DblNumMat Merr(Mext.m(), Mext.n());
    for(int a=0; a<Mext.m(); a++)
      for(int b=0; b<Mext.n(); b++)
        Merr(a,b) = Mext(a,b) - Mapp(a,b);
    cerr<<"rel err "<<sqrt(energy(Merr))/sqrt(energy(Mext))<<endl;
  }
  return 0;
}

//complex version
int ddlowrank(int m, int n, int (*sample)(vector<int>&, vector<int>&, ZpxNumMat&), double eps, int npk, 
	    vector<int>& cidx, vector<int>& ridx, ZpxNumMat& mid)
{
  iA(m>0 && n>0);
  {
    int nc = min(npk,n);
    vector<int> cs(nc);
    for(int k=0; k<nc; k++)      cs[k] = int( floor(drand48()*(n-nc)) );
    sort(cs.begin(), cs.end());
    for(int k=0; k<nc; k++)      cs[k] += k;
    vector<int> rs(m);
    for(int k=0; k<m; k++)      rs[k] = k;
    //
    ZpxNumMat M2tmp;    iC( (*sample)(rs, cs, M2tmp) );
    ZpxNumMat M2(M2tmp.n(), M2tmp.m());
    for(int i=0; i<M2tmp.m(); i++)      for(int j=0; j<M2tmp.n(); j++)	M2(j,i) = conj(M2tmp(i,j));
    //
    int m = M2.m();
    int n = M2.n();
    int lda = m;
    NumVec<int> jpvt(n);    setvalue(jpvt, int(0));
    ZpxNumVec tau(max(m,n));
    ZpxNumVec work(3*n);
    DblNumVec rwork(6*n);
    int lwork = 6*n;
    int info;
    zgeqp3_(&m, &n, 
	    (MKL_Complex16*) M2.data(), &lda, jpvt.data(), 
	    (MKL_Complex16*) tau.data(), 
	    (MKL_Complex16*) work.data(), &lwork, rwork.data(), &info);    
    iA(info==0);
    double cutoff = eps*abs(M2(0,0));
    int cnt=0;
    for(int k=0; k<min(m,n); k++)      if(abs(M2(k,k))>cutoff)	cnt++;
    ridx.resize(cnt);
    for(int k=0; k<cnt; k++)      ridx[k] = (jpvt(k)-1);
    cerr<<"ROWS "; for(int k=0; k<cnt; k++)      cerr<<ridx[k]<<" ";    cerr<<endl;
  }
  {
    int nr = min(npk,m);
    vector<int> rs(nr);
    for(int k=0; k<nr; k++)      rs[k] = int( floor(drand48()*(m-nr)) );
    sort(rs.begin(), rs.end());
    for(int k=0; k<nr; k++)      rs[k] += k;
    //
    for(unsigned int k=0; k<ridx.size(); k++)      rs.push_back(ridx[k]);
    sort(rs.begin(), rs.end());
    vector<int>::iterator newend = unique(rs.begin(), rs.end());
    rs.resize(newend-rs.begin());
    //
    vector<int> cs(n);
    for(int k=0; k<n; k++)      cs[k] = k;
    //
    ZpxNumMat M1;    iC( (*sample)(rs, cs, M1) );
    //
    int m = M1.m();
    int n = M1.n();
    int lda = m;
    NumVec<int> jpvt(n);      setvalue(jpvt, int(0));
    ZpxNumVec tau(max(m,n));
    ZpxNumVec work(3*n);
    int lwork = 6*n;
    DblNumVec rwork(6*n);
    int info;
    zgeqp3_(&m, &n, 
	    (MKL_Complex16*) M1.data(), &lda, jpvt.data(), 
	    (MKL_Complex16*) tau.data(), 
	    (MKL_Complex16*) work.data(), &lwork, rwork.data(), &info);    
    iA(info==0);
    double cutoff = eps*abs(M1(0,0)); //the diagonal element
    int cnt=0;
    for(int k=0; k<min(m,n); k++)	if(abs(M1(k,k))>cutoff)	  cnt++;
    cidx.resize(cnt);
    for(int k=0; k<cnt; k++)	cidx[k] = (jpvt(k)-1);
    cerr<<"COLS "; for(int k=0; k<cnt; k++)      cerr<<cidx[k]<<" ";    cerr<<endl;
  }
  {
    int nc = min(npk,n);
    vector<int> cs(nc);
    for(int k=0; k<nc; k++)      cs[k] = int( floor(drand48()*(n-nc)) );
    sort(cs.begin(), cs.end());
    for(int k=0; k<nc; k++)      cs[k] += k;
    for(unsigned int k=0; k<cidx.size(); k++)      cs.push_back(cidx[k]);
    sort(cs.begin(), cs.end());
    vector<int>::iterator csnewend = unique(cs.begin(), cs.end());
    cs.resize(csnewend-cs.begin());
    //
    int nr = min(npk,m);
    vector<int> rs(nr);
    for(int k=0; k<nr; k++)      rs[k] = int( floor(drand48()*(m-nr)) );
    sort(rs.begin(), rs.end());
    for(int k=0; k<nr; k++)      rs[k] += k;
    for(unsigned int k=0; k<ridx.size(); k++)      rs.push_back(ridx[k]);
    sort(rs.begin(), rs.end());
    vector<int>::iterator rsnewend = unique(rs.begin(), rs.end());
    rs.resize(rsnewend-rs.begin());
    //
    ZpxNumMat M1;    iC( (*sample)(rs,cidx,M1) );
    ZpxNumMat IM1;    iC( ddpinv(M1, 1e-16, IM1) );
    ZpxNumMat M2;    iC( (*sample)(ridx,cs,M2) );
    ZpxNumMat IM2;    iC( ddpinv(M2, 1e-16, IM2) );
    ZpxNumMat M3;    iC( (*sample)(rs,cs,M3) );
    ZpxNumMat tmp(M3.m(), IM2.n());
    iC( zzgemm(1.0, M3, IM2, 0.0, tmp) );
    mid.resize(IM1.m(), tmp.n());
    iC( zzgemm(1.0, IM1, tmp, 0.0, mid) );
  }
  if(1) {
    int nc = min(npk,n);
    vector<int> cs(nc);
    for(int k=0; k<nc; k++)      cs[k] = int( floor(drand48()*n) );
    int nr = min(npk,m);
    vector<int> rs(nr);
    for(int k=0; k<nr; k++)      rs[k] = int( floor(drand48()*m) );
    ZpxNumMat M1;
    iC( (*sample)(rs,cidx,M1) );
    ZpxNumMat M2;
    iC( (*sample)(ridx,cs,M2) );
    ZpxNumMat Mext;
    iC( (*sample)(rs,cs,Mext) );
    ZpxNumMat Mapp(rs.size(), cs.size());
    ZpxNumMat tmp(mid.m(),M2.n());
    iC( zzgemm(1.0, mid, M2, 0.0, tmp) );
    iC( zzgemm(1.0, M1, tmp, 0.0, Mapp) );
    ZpxNumMat Merr(Mext.m(), Mext.n());
    for(int a=0; a<Mext.m(); a++)
      for(int b=0; b<Mext.n(); b++)
	Merr(a,b) = Mext(a,b) - Mapp(a,b);
    cerr<<"lowrank rel err "<<sqrt(energy(Merr))/sqrt(energy(Mext))<<endl;
  }
  return 0;
}













