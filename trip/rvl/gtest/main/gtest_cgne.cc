/*************************************************************************

Copyright Rice University, 2011.
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/
#include "gtest/gtest.h"
#include "rnmat.hh"
#include "adjtest.hh"
#include "cgnealg.hh"

#define RANDSEED 19490615

namespace {

  using namespace RVL;
  using namespace RVLAlg;
  using namespace RVLUmin;

  template<typename T>
  int CGNEMain(int rows, 
	       int cols, 
	       typename ScalarFieldTraits<T>::AbsType _nlvl, 
	       int maxit, 
	       typename ScalarFieldTraits<T>::AbsType maxsv,
	       typename ScalarFieldTraits<T>::AbsType opnoise,
	       typename ScalarFieldTraits<T>::AbsType rtol, 
	       typename ScalarFieldTraits<T>::AbsType nrtol, 
	       typename ScalarFieldTraits<T>::AbsType Delta, 
	       bool precond,
	       string output, 
	       ostream & str) {
    
    typedef typename ScalarFieldTraits<T>::AbsType atype;

    try {
      // build output stream
      ofstream strout(output.c_str());
      strout<<scientific;

      //      str<<"1. Construct domain, range spaces"<<endl;
      RnSpace<T> dom(cols);
      RnSpace<T> rng(rows);
      //      str<<endl;

      //      str<<"2. Build LinearOp (GenMat), do adjoint test\n";
      GenMat<T> A(dom,rng);
      GenMat<T> M(dom,dom);
#ifdef RANDSEED
      int randseed=RANDSEED;
#else
      int randseed=getpid();
#endif
      RVLRandomize<T> rndop(randseed,-opnoise,opnoise);
      A.eval(rndop);
      RVLAssignConst<T> ac(ScalarFieldTraits<T>::Zero());
      M.eval(ac);
      for (int i=0;i<min(cols,rows);i++) {
	A.getElement(i,i)+=ScalarFieldTraits<atype>::One() + (i*(maxsv-ScalarFieldTraits<atype>::One()))/min(cols,rows);
	M.getElement(i,i)=ScalarFieldTraits<atype>::One()/(ScalarFieldTraits<atype>::One() + (i*(maxsv-ScalarFieldTraits<atype>::One()))/min(cols,rows));
	M.getElement(i,i) *= M.getElement(i,i);
	//	cerr<<"i="<<i<<" A[i,i]="<<A.getElement(i,i)<<endl;
      }
      /*
      if (!AdjointTest<T>(A,rnd,strout)) {
	str<<"*** adjoint test failed - result in "<<output<<"; exit\n";
	return INT_MAX;
      }
      else {
	str<<"*** adjoint test passed - result in "<<output<<"; proceed\n";
      }
      str<<endl;
      */

      //      str<<"3. Build sol, rhs vectors\n";
      RVLRandomize<T> rnd(3*randseed,-1.0,1.0);
      Vector<T> sol(dom);
      Vector<T> est(dom);
      // initialize exact solution to random vector
      sol.eval(rnd);
      // initialize estimated solution to zero
      est.zero();
      // create rhs, initialize to noisy version of exact sol image
      Vector<T> rhs(rng);
      A.applyOp(sol,rhs);
      atype rhsnorm = rhs.norm();
      Vector<T> err(rng);
      err.eval(rnd);
      T nlvl = _nlvl*rhs.norm()/err.norm();
      rhs.linComb(nlvl,err);
      //      str<<endl;

      //      str<<"4. Build, run CGNE algorithm \n";
      //      str<<endl;    
      strout<<endl;
      strout<<"********************************************************"<<endl;
      strout<<"      Conjugate Gradient Algorithm for Normal Eqns"<<endl;
      strout<<endl;
      strout<<"max iterations               = "<<maxit<<endl;
      strout<<"residual tolerance           = "<<rtol<<endl;
      strout<<"normal res tolerance         = "<<nrtol<<endl;
      strout<<"trust radius                 = "<<Delta<<endl;
      strout<<"singular value range         = [1, "<<maxsv<<"]"<<endl;
      strout<<"matrix noise level           = "<<opnoise<<endl;
      strout<<"data relerr                  = "<<_nlvl<<endl;
      if (precond) 
	strout<<"preconditioned with inverse diagonal squared\n";
      strout<<endl;
      strout<<"**************** iteration summary *********************"<<endl;

      // allocate rnorm, nrnorm scalars
      atype rnorm = numeric_limits<atype>::max();
      atype nrnorm = numeric_limits<atype>::max();
      // construct CGNEAlg
      int actit=0;
      if (precond) {
        CGNEAlg<T> iter(est,A,M,rhs,rnorm,nrnorm,rtol,nrtol,maxit,Delta,strout);
	iter.run();
	actit=iter.getCount();
      }
      else {
	CGNEAlg<T> iter(est,A,rhs,rnorm,nrnorm,rtol,nrtol,maxit,Delta,strout);	
	iter.run();
	actit=iter.getCount();
      }

      // display results
      atype snorm = sol.norm();
      atype enorm = est.norm();
      sol.linComb(-ScalarFieldTraits<atype>::One(),est);
      atype relerr=ScalarFieldTraits<atype>::One();
      if (ProtectedDivision<atype>(sol.norm(),snorm,relerr)) {
	RVLException e;
	e<<"a-oo-ga! a-oo-ga! solution was zero!!!\n";
	throw e;
      }
      strout<<"actual number of iterations  = "<<actit;
      strout<<"norm of noise-free data      = "<<rhsnorm<<endl;
      strout<<"norm of noisy data           = "<<rhs.norm()<<endl;
      strout<<"residual norm                = "<<rnorm<<endl;
      strout<<"residual relerr              = "<<rnorm/rhsnorm<<endl;
      strout<<"gradient norm                = "<<nrnorm<<endl;
      strout<<"norm of solution             = "<<enorm<<endl;
      strout<<"norm of target vector        = "<<snorm<<endl;
      strout<<"relative error in solution   = "<<relerr<<endl;
      strout.close();

      return actit;
    }
    catch (RVLException & e) {
      e.write(str);
      exit(1);
    }
  }

  class CGNETest: public ::testing::Test {
    
  public:
    CGNETest() {}
  };

  TEST_F(CGNETest, real_flt_mn0_dn0) {
    typedef float T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0f;
    int maxit=100; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0f;
    ScalarFieldTraits<T>::AbsType opnoise=0.0f;
    ScalarFieldTraits<T>::AbsType rtol=0.01f;
    ScalarFieldTraits<T>::AbsType nrtol=0.01f;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0f;
    bool precond=false;
    string output="report_real_flt_mn0_dn0.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, real_dbl_mn0_dn0) {
    typedef double T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0;
    int maxit=100; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0;
    ScalarFieldTraits<T>::AbsType opnoise=0.0;
    ScalarFieldTraits<T>::AbsType rtol=0.01;
    ScalarFieldTraits<T>::AbsType nrtol=0.01;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0;
    bool precond=false;
    string output="report_real_dbl_mn0_dn0.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, complex_flt_mn0_dn0) {
    typedef complex<float> T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0;
    int maxit=100; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0;
    ScalarFieldTraits<T>::AbsType opnoise=0.0;
    ScalarFieldTraits<T>::AbsType rtol=0.01;
    ScalarFieldTraits<T>::AbsType nrtol=0.01;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0;
    bool precond=false;
    string output="report_complex_flt_mn0_dn0.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, complex_dbl_mn0_dn0) {
    typedef complex<double> T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0;
    int maxit=100; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0;
    ScalarFieldTraits<T>::AbsType opnoise=0.0;
    ScalarFieldTraits<T>::AbsType rtol=0.01;
    ScalarFieldTraits<T>::AbsType nrtol=0.01;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0;
    bool precond=false;
    string output="report_complex_dbl_mn0_dn0.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, real_dbl_mn0_dn0_tr) {
    typedef double T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0;
    int maxit=100; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0;
    ScalarFieldTraits<T>::AbsType opnoise=0.0;
    ScalarFieldTraits<T>::AbsType rtol=0.01;
    ScalarFieldTraits<T>::AbsType nrtol=0.01;
    ScalarFieldTraits<T>::AbsType Delta=1.0;
    bool precond=false;
    string output="report_real_dbl_mn0_dn0_tr.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, real_dbl_mn0_dn0_2_tr) {
    typedef double T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.2;
    int maxit=100; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0;
    ScalarFieldTraits<T>::AbsType opnoise=0.0;
    ScalarFieldTraits<T>::AbsType rtol=0.01;
    ScalarFieldTraits<T>::AbsType nrtol=0.01;
    ScalarFieldTraits<T>::AbsType Delta=1.0;
    bool precond=false;
    string output="report_real_dbl_mn0_dn0.2_tr.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, real_flt_mn0_dn0_pre) {
    typedef float T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0f;
    int maxit=5; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0f;
    ScalarFieldTraits<T>::AbsType opnoise=0.0f;
    ScalarFieldTraits<T>::AbsType rtol=0.01f;
    ScalarFieldTraits<T>::AbsType nrtol=0.01f;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0f;
    bool precond=true;
    string output="report_real_flt_mn0_dn0_pre.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, real_dbl_mn0_dn0_pre) {
    typedef double T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0;
    int maxit=5; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0;
    ScalarFieldTraits<T>::AbsType opnoise=0.0;
    ScalarFieldTraits<T>::AbsType rtol=0.01;
    ScalarFieldTraits<T>::AbsType nrtol=0.01;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0;
    bool precond=true;
    string output="report_real_dbl_mn0_dn0_pre.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, complex_flt_mn0_dn0_pre) {
    typedef complex<float> T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0;
    int maxit=5; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0;
    ScalarFieldTraits<T>::AbsType opnoise=0.0;
    ScalarFieldTraits<T>::AbsType rtol=0.01;
    ScalarFieldTraits<T>::AbsType nrtol=0.01;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0;
    bool precond=true;
    string output="report_complex_flt_mn0_dn0_pre.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, complex_dbl_mn0_dn0_pre) {
    typedef complex<double> T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0;
    int maxit=5; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0;
    ScalarFieldTraits<T>::AbsType opnoise=0.0;
    ScalarFieldTraits<T>::AbsType rtol=0.01;
    ScalarFieldTraits<T>::AbsType nrtol=0.01;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0;
    bool precond=true;
    string output="report_complex_dbl_mn0_dn0_pre.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, real_flt_mnem3_dn0_pre) {
    typedef float T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0f;
    int maxit=5; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0f;
    ScalarFieldTraits<T>::AbsType opnoise=0.001f;
    ScalarFieldTraits<T>::AbsType rtol=0.01f;
    ScalarFieldTraits<T>::AbsType nrtol=0.01f;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0f;
    bool precond=true;
    string output="report_real_flt_mnem3_dn0_pre.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, real_dbl_mnem3_dn0_pre) {
    typedef double T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0;
    int maxit=5; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0;
    ScalarFieldTraits<T>::AbsType opnoise=0.001;
    ScalarFieldTraits<T>::AbsType rtol=0.01;
    ScalarFieldTraits<T>::AbsType nrtol=0.01;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0;
    bool precond=true;
    string output="report_real_dbl_mnem3_dn0_pre.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, complex_flt_mnem3_dn0_pre) {
    typedef complex<float> T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0;
    int maxit=5; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0;
    ScalarFieldTraits<T>::AbsType opnoise=0.001;
    ScalarFieldTraits<T>::AbsType rtol=0.01;
    ScalarFieldTraits<T>::AbsType nrtol=0.01;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0;
    bool precond=true;
    string output="report_complex_flt_mnem3_dn0_pre.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, complex_dbl_mnem3_dn0_pre) {
    typedef complex<double> T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0;
    int maxit=5; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0;
    ScalarFieldTraits<T>::AbsType opnoise=0.001;
    ScalarFieldTraits<T>::AbsType rtol=0.01;
    ScalarFieldTraits<T>::AbsType nrtol=0.01;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0;
    bool precond=true;
    string output="report_complex_dbl_mnem3_dn0_pre.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, real_flt_mnem2_dn0_pre) {
    typedef float T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0f;
    int maxit=5; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0f;
    ScalarFieldTraits<T>::AbsType opnoise=0.01f;
    ScalarFieldTraits<T>::AbsType rtol=0.01f;
    ScalarFieldTraits<T>::AbsType nrtol=0.01f;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0f;
    bool precond=true;
    string output="report_real_flt_mnem2_dn0_pre.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, real_dbl_mnem2_dn0_pre) {
    typedef double T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0;
    int maxit=5; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0;
    ScalarFieldTraits<T>::AbsType opnoise=0.01;
    ScalarFieldTraits<T>::AbsType rtol=0.01;
    ScalarFieldTraits<T>::AbsType nrtol=0.01;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0;
    bool precond=true;
    string output="report_real_dbl_mnem2_dn0_pre.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, complex_flt_mnem2_dn0_pre) {
    typedef complex<float> T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0;
    int maxit=8; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0;
    ScalarFieldTraits<T>::AbsType opnoise=0.01;
    ScalarFieldTraits<T>::AbsType rtol=0.01;
    ScalarFieldTraits<T>::AbsType nrtol=0.01;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0;
    bool precond=true;
    string output="report_complex_flt_mnem2_dn0_pre.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, complex_dbl_mnem2_dn0_pre) {
    typedef complex<double> T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0;
    int maxit=8; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0;
    ScalarFieldTraits<T>::AbsType opnoise=0.01;
    ScalarFieldTraits<T>::AbsType rtol=0.01;
    ScalarFieldTraits<T>::AbsType nrtol=0.01;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0;
    bool precond=true;
    string output="report_complex_dbl_mnem2_dn0_pre.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, real_flt_mnem1_dn0_pre) {
    typedef float T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0f;
    int maxit=10; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0f;
    ScalarFieldTraits<T>::AbsType opnoise=0.1f;
    ScalarFieldTraits<T>::AbsType rtol=0.01f;
    ScalarFieldTraits<T>::AbsType nrtol=0.01f;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0f;
    bool precond=true;
    string output="report_real_flt_mnem1_dn0_pre.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, real_dbl_mnem1_dn0_pre) {
    typedef double T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0;
    int maxit=10; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0;
    ScalarFieldTraits<T>::AbsType opnoise=0.1;
    ScalarFieldTraits<T>::AbsType rtol=0.01;
    ScalarFieldTraits<T>::AbsType nrtol=0.01;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0;
    bool precond=true;
    string output="report_real_dbl_mnem1_dn0_pre.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, complex_flt_mnem1_dn0_pre) {
    typedef complex<float> T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0;
    int maxit=15; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0;
    ScalarFieldTraits<T>::AbsType opnoise=0.1;
    ScalarFieldTraits<T>::AbsType rtol=0.01;
    ScalarFieldTraits<T>::AbsType nrtol=0.01;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0;
    bool precond=true;
    string output="report_complex_flt_mnem1_dn0_pre.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

  TEST_F(CGNETest, complex_dbl_mnem1_dn0_pre) {
    typedef complex<double> T;
    int rows=261;
    int cols=121;
    ScalarFieldTraits<T>::AbsType _nlvl=0.0;
    int maxit=15; 
    ScalarFieldTraits<T>::AbsType maxsv=10.0;
    ScalarFieldTraits<T>::AbsType opnoise=0.1;
    ScalarFieldTraits<T>::AbsType rtol=0.01;
    ScalarFieldTraits<T>::AbsType nrtol=0.01;
    ScalarFieldTraits<T>::AbsType Delta=1000000.0;
    bool precond=true;
    string output="report_complex_dbl_mnem1_dn0_pre.txt"; 
    ostream & str = cerr;
    int actit = CGNEMain<T>(rows,cols,_nlvl,maxit,maxsv,opnoise,rtol,nrtol,Delta,precond,output,str);
    EXPECT_LT(actit,maxit);
  }

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  int err = RUN_ALL_TESTS();
  return err;
}
