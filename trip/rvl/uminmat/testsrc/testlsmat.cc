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


#include "rnmat.hh"
#include "adjtest.hh"
#include "TRLSMat.hh"

using namespace RVL;
using namespace RVLAlg;
using namespace RVLUmin;

/** Umin Unit Test 1
    Purpose: basic TRLS test - scalar type = float
*/

template<typename T>
void LSMatTest(int rows, int cols,  
	       typename ScalarFieldTraits<T>::AbsType _nlvl, 
	       typename ScalarFieldTraits<T>::AbsType Delta, 
	       string testno, string datatype, string output, ostream & str) {

  typedef typename ScalarFieldTraits<T>::AbsType atype;

  try {
    str<<endl;
    str<<"/****************************************************"<<endl;
    str<<" *            BEGIN UMin UNIT TEST "<<testno<<endl;
    str<<" * Purpose: test LAPACK-based LS solver"<<endl;
    str<<" *   scalar type    = "<<datatype<<endl;
    str<<" *   noise level    = "<<_nlvl<<endl;
    str<<" *   trust radius   = "<<Delta<<endl;
    str<<" *   results reported in "<<output<<endl;
    str<<" ****************************************************/"<<endl; 
    str<<endl;
    ofstream strout(output.c_str());
    strout<<scientific;

    str<<"1. Construct domain, range spaces"<<endl;
    RnSpace<T> dom(cols);
    RnSpace<T> rng(rows);
    str<<endl;
    str<<"2. Build GenMat\n";
    GenMat<T> A(dom,rng);
    RVLRandomize<T> rnd;
    A.eval(rnd);
    for (int i=0;i<min(cols,rows);i++) 
      A.getElement(i,i)+=ScalarFieldTraits<T>::One();
    str<<endl;
    str<<"3. Build sol, rhs\n";
    Vector<T> sol(dom);
    Vector<T> est(dom);
    sol.eval(rnd);
    est.zero();
    Vector<T> rhs(rng);
    A.applyOp(sol,rhs);
    atype rhsnorm = rhs.norm();
    Vector<T> err(rng);
    err.eval(rnd);
    // assign atype to Scalar, double to atype - should work in all cases
    T nlvl = _nlvl*rhs.norm()/err.norm();
    rhs.linComb(nlvl,err);
    str<<endl;
    str<<"4. Build, run TRLSMat algorithm \n";
    str<<endl;    
    atype rnorm = numeric_limits<atype>::max();
    atype nrnorm = numeric_limits<atype>::max();
    str<<endl;    
    strout<<endl<<"*************************************************************"<<endl;
    strout<<"* Test of LAPACK-based least squares solver with TR truncation"<<endl;
    strout<<"*   noise level          = "<<_nlvl<<endl;
    strout<<"*   trust radius         = "<<Delta<<endl;
    strout<<"************************************************************"<<endl;

    LSMatAlg<T> iter(est,A,rhs,rnorm,nrnorm,Delta,strout);
    iter.run();
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

int main() {

  try {
    PlantSeeds(19490615);
    int rows=261;
    int cols=101;

    LSMatTest<float>(rows, cols, 0.0, 1000000.0,"3.1","float","lsflt.rpt",cout);  
    LSMatTest<double>(rows, cols, 0.1, 1000000.0,"3.2","double","lsdbl.rpt",cout);  
    LSMatTest<complex<float> >(rows, cols, 0.1, 1000000.0,"3.3","complex<float>","lsfltcplx.rpt",cout);  
    LSMatTest<complex<double> >(rows, cols, 0.5, 1000000.0,"3.4","complex<double>","lsdblcplx.rpt",cout);  
    LSMatTest<float>(rows, cols, 0.0, 1.0,"3.5","float","lsflt1.rpt",cout);  
    LSMatTest<double>(rows, cols, 0.1, 1.0,"3.6","double","lsdbl1.rpt",cout);  
  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
