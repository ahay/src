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

namespace {


  using namespace RVL;
  using namespace RVLAlg;
  using namespace RVLUmin;

  /** Umin Unit Test 1
      Purpose: basic CG alg - scalar type = float
  */

  template<typename T>
  void CGNEMain(int rows, int cols, typename ScalarFieldTraits<T>::AbsType _nlvl, 
		int maxit, 
		typename ScalarFieldTraits<T>::AbsType rtol, 
		typename ScalarFieldTraits<T>::AbsType nrtol, 
		typename ScalarFieldTraits<T>::AbsType Delta, 
		string testno, string datatype, string output, 
		ostream & str) {

    typedef typename ScalarFieldTraits<T>::AbsType atype;

    try {
      str<<endl;
      str<<"/****************************************************"<<endl;
      str<<" *            BEGIN UMin UNIT TEST "<<testno<<endl;
      str<<" * Purpose: test CGNE alg scalar type ="<<datatype<<endl;
      str<<" *   noise level = "<<_nlvl<<endl;
      str<<" *   results reported in "<<output<<endl;
      str<<" ****************************************************/"<<endl; 
      str<<endl;

      // build output stream
      ofstream strout(output.c_str());
      strout<<scientific;

      str<<"1. Construct domain, range spaces"<<endl;
      RnSpace<T> dom(cols);
      RnSpace<T> rng(rows);
      str<<endl;

      str<<"2. Build LinearOp (GenMat), do adjoint test\n";
      GenMat<T> A(dom,rng);
      RVLRandomize<T> rnd;
      A.eval(rnd);
      for (int i=0;i<min(cols,rows);i++) 
	A.getElement(i,i)+=ScalarFieldTraits<T>::One();
      if (AdjointTest<T>(A,rnd,strout)) {
	str<<"*** adjoint test passed - result in "<<output<<"; proceed\n";
      }
      else {
	str<<"*** adjoint test failed - result in "<<output<<"; exit\n";
	return;
      }
      str<<endl;

      str<<"3. Build sol, rhs vectors\n";
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
      str<<endl;

      str<<"4. Build, run CGNE algorithm \n";
      str<<endl;    
      strout<<endl<<"*******************************************************"<<endl;
      strout<<"* Test of Conjugate Gradient Algorithm for Normal Eqns"<<endl;
      strout<<"* max iterations       = "<<maxit<<endl;
      strout<<"* residual tolerance   = "<<rtol<<endl;
      strout<<"* normal res tolerance = "<<nrtol<<endl;
      strout<<"* trust radius         = "<<Delta<<endl;
      strout<<"*******************************************************"<<endl;

      // allocate rnorm, nrnorm scalars
      atype rnorm = numeric_limits<atype>::max();
      atype nrnorm = numeric_limits<atype>::max();
      // construct CGNEAlg
      CGNEAlg<T> iter(est,A,rhs,rnorm,nrnorm,rtol,nrtol,maxit,Delta,strout);
      // run it
      iter.run();

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
      strout<<"\n ******* summary ********  "<<endl;
      strout<<"norm of noise-free data    = "<<rhsnorm<<endl;
      strout<<"norm of noisy data         = "<<rhs.norm()<<endl;
      strout<<"residual norm              = "<<rnorm<<endl;
      strout<<"residual relerr            = "<<rnorm/rhsnorm<<endl;
      strout<<"data relerr                = "<<_nlvl<<endl;
      strout<<"gradient norm              = "<<nrnorm<<endl;
      strout<<"norm of solution           = "<<enorm<<endl;
      strout<<"norm of target vector      = "<<snorm<<endl;
      strout<<"relative error in solution = "<<relerr<<endl;
      strout.close();
    }
    catch (RVLException & e) {
      e.write(str);
      exit(1);
    }
  }

  class CGNETest: public ::testing::Test {
    
  public:
   
    int rows;
    int cols;
    ofstream dump;
    CGNETest(): rows(261), cols(121), dump("dump.txt") {
      PlantSeeds(19490615);
    }
  };

  TEST_F(CGNETest, real_flt) {
    CGNEMain<float>(rows,cols,0.0f,100,0.001f,0.001f,1000000.0f,"2.1","float","./cgneflt.rpt",dump);
  }
  TEST_F(CGNETest, real_dbl) {
    CGNEMain<double>(rows,cols,0.1,100,0.001,0.001,1000000.0,"2.2","double","./cgnedbl.rpt",dump);
  }
  TEST_F(CGNETest, complex_float) {
    CGNEMain<complex<float> >(rows,cols,0.1,100,0.001,0.001,1000000.0,"2.3","complex<float>","./cgnefltcplx.rpt",dump);
  }
  TEST_F(CGNETest, complex_double) {
    CGNEMain<double>(rows,cols,0.2,100,0.001,0.001,1000000.0,"2.4","complex<double>","./cgnedblcplx.rpt",dump);
  }
  TEST_F(CGNETest, real_double_trustrad) {
    CGNEMain<float>(rows,cols,0.0,100,0.001,0.001,1.0,"2.5","float","./cgneflt1.rpt",dump);    
  }
  TEST_F(CGNETest, real_double_noise_trustrad) {
    CGNEMain<double>(rows,cols,0.2,100,0.001,0.001,1.0,"2.6","complex<double>","./cgnedblcplx1.rpt",dump);
  }

}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  int err = RUN_ALL_TESTS();
  return err;
}
