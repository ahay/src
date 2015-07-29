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
#include "chebalg.hh"

using namespace RVL;
using namespace RVLAlg;
using namespace RVLUmin;

/** Umin Unit Test 1
    Purpose: basic CG alg - scalar type = float
*/

template<typename T>
void ChebTest(int rows, 
	      int cols, 
	      typename ScalarFieldTraits<T>::AbsType _nlvl, // noise level
	      int maxit, 
	      typename ScalarFieldTraits<T>::AbsType minsig, // max singular val
	      typename ScalarFieldTraits<T>::AbsType maxsig, // max singular val
	      typename ScalarFieldTraits<T>::AbsType gamma,  // inversion level
              typename ScalarFieldTraits<T>::AbsType epsilon,// error reduction
              typename ScalarFieldTraits<T>::AbsType alpha,  // 'fudge factor'
	      string testno, string datatype, string output, 
	      ostream & str) {

  typedef typename ScalarFieldTraits<T>::AbsType atype;

  try {
    str<<endl;
    str<<"/****************************************************"<<endl;
    str<<" *            BEGIN UMin UNIT TEST "<<testno<<endl;
    str<<" * Purpose: test Cheb alg scalar type ="<<datatype<<endl;
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
    int n = min(cols,rows);
    atype one = ScalarFieldTraits<atype>::One();
    for (int i=0;i<n;i++) {
      atype t=((atype)i)/((atype)n);
      A.getElement(i,i)=sqrt((one-t)*minsig*minsig + t*maxsig*maxsig);
    }
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
    Vector<T> nrhs(rng);
    A.applyAdjOp(rhs,nrhs);
    atype nrhsnorm = nrhs.norm();
    str<<endl;

    str<<"4. Build, run Cheb algorithm \n";
    str<<endl;    

    // allocate rnorm, nrnorm scalars
    atype rnorm = numeric_limits<atype>::max();
    atype nrnorm = numeric_limits<atype>::max();
    ChebAlg<T> iter(est,A,rhs,rnorm,nrnorm,gamma,epsilon,alpha,maxit,strout);
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
      
    Vector<T> ATb(A.getDomain());
    A.applyAdjOp(rhs,ATb);
    rhsnorm = ATb.norm();
    strout<<"*******************************************************"<<endl;
    strout<<"*     Test of Chebyshev iteration for Normal Eqns     *"<<endl;
    strout<<endl<<"****** Example Attributes  "<<endl;
    strout<<"max iterations               = "<<maxit<<endl;
    strout<<"data noise level             = "<<_nlvl<<endl;
    strout<<"min singular val             = "<<minsig<<endl;
    strout<<"max singular val             = "<<maxsig<<endl;
    strout<<"rcond for normal eqns        = "<<minsig*minsig/(maxsig*maxsig)<<endl;
    strout<<"inversion level              = "<<gamma<<endl;
    strout<<"target normal res redn       = "<<epsilon<<endl;
    strout<<"fudge factor                 = "<<alpha<<endl;

    strout<<endl<<"****** Iteration Summary "<<endl;
    strout<<"iterations since restart     = "<<iter.getCount()<<endl;
    strout<<"total iterations             = "<<iter.getTotalCount()<<endl;
    strout<<"number of restarts           = "<<iter.getRestartCount()<<endl;
    strout<<"norm of ATb                  = "<<nrhsnorm<<endl;
    strout<<"normal res norm              = "<<nrnorm<<endl;
    strout<<"actual normal res redn       = "<<nrnorm/nrhsnorm<<endl;
    strout<<"data relerr                  = "<<_nlvl<<endl;
    strout<<"data norm                    = "<<rhsnorm<<endl;
    strout<<"res norm                     = "<<rnorm<<endl;
    strout<<"actual res redn              = "<<rnorm/rhsnorm<<endl;
    strout<<"norm of exact solution       = "<<snorm<<endl;
    strout<<"norm of estimated solution   = "<<enorm<<endl;
    strout<<"relative error in solution   = "<<relerr<<endl;
    strout<<"final value of LS fcn        = ";
    Vector<T> val(A.getRange());
    A.applyOp(est,val);
    val.linComb(-1.0,rhs);
    strout<< 1.0/2.0*val.normsq()<<endl;
    strout<<"*******************************************************"<<endl;
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

    int cols=100;

    ChebTest<double>(cols,cols,0.0,1000,0.1,1.0,0.04,0.001,1.1,
		     "0.0","double","./Chebdbl0.rpt",cout);
    ChebTest<double>(cols,cols,0.0,1000,0.1,1.0,0.01,0.001,1.1,
		     "0.1","double","./Chebdbl1.rpt",cout);
    ChebTest<double>(cols,cols,0.0,1000,0.01,1.0,0.04,0.001,1.1,
		     "0.2","double","./Chebdbl2.rpt",cout);
    ChebTest<double>(cols,cols,0.0,1000,0.01,1.0,0.0001,0.001,1.1,
		     "0.3","double","./Chebdbl3.rpt",cout);
    ChebTest<double>(cols,cols,0.0,1000,0.0,1.0,0.04,0.001,1.1,
		     "0.4","double","./Chebdbl4.rpt",cout);
    ChebTest<double>(cols,cols,0.0,1000,0.0,1.0,0.0001,0.001,1.1,
		     "0.5","double","./Chebdbl5.rpt",cout);
  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
