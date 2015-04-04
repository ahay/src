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
#include "cgalg.hh"

using namespace RVL;
using namespace RVLAlg;
using namespace RVLUmin;

/** Umin Unit Test 1
    Purpose: basic CG alg - scalar type = float
*/

template<typename T>
void CGTest(int rows, int cols, 
	    int maxit, typename ScalarFieldTraits<T>::AbsType rtol, 
	    string testno, string datatype, string output, 
	    ostream & str) {

  typedef typename ScalarFieldTraits<T>::AbsType atype;

  try {
    str<<endl;
    str<<"/****************************************************"<<endl;
    str<<" *            BEGIN UMin UNIT TEST "<<testno<<endl;
    str<<" * Purpose: test basic cg alg scalar type ="<<datatype<<endl;
    str<<" *   results reported in "<<output<<endl;
    str<<" ****************************************************/"<<endl; 
    str<<endl;
    ofstream strout(output.c_str());
    strout<<scientific;
    str<<"1. Construct domain, range spaces"<<endl;
    RnSpace<T> dom(cols);
    RnSpace<T> rng(rows);
    str<<endl;
    str<<"2. Build GenMat, do adjoint test\n";
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
    str<<endl<<"3. Build Normal operator - will be 121x121\n";
    NormalLinearOp<T> N(A);
    str<<endl;
    str<<"4. Build sol, rhs\n";
    Vector<T> sol(dom);
    Vector<T> est(dom);
    sol.eval(rnd);
    est.zero();
    Vector<T> rhs(dom);
    N.applyOp(sol,rhs);
    str<<endl;
    str<<"5. Build, run CG algorithm - results reported in "<<output<<"\n";
    str<<endl;    
    strout<<endl<<"*************************************************"<<endl;
    strout<<"* Test of Conjugate Gradient Algorithm"<<endl;
    strout<<"* max iterations      = "<<maxit<<endl;
    strout<<"* residual tolerance  = "<<rtol<<endl;
    strout<<"*************************************************"<<endl;
    atype rnormsq = numeric_limits<atype>::max();
    atype maxstep = numeric_limits<atype>::max();
    CGAlg<T> iter(est,N,rhs,rnormsq,rtol,maxit,maxstep,strout);
    iter.run();
    str<<"converged in "<<iter.getCount()<<" steps"<<endl;
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
    strout<<"norm of data                 = "<<rhs.norm()<<endl;
    strout<<"residual norm                = "<<sqrt(rnormsq)<<endl;
    strout<<"residual relerr              = "<<sqrt(rnormsq)/rhs.norm()<<endl;
    strout<<"norm of solution             = "<<enorm<<endl;
    strout<<"norm of target vector        = "<<snorm<<endl;
    strout<<"relative error in solution   = "<<relerr<<endl;
    
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
    int cols=121;

    CGTest<float>(rows,cols,100, 0.001,"1.1","float","./cgflt.txt",cout);
    CGTest<double>(rows,cols,100, 0.001,"1.2","double","./cgdbl.txt",cout);
    CGTest< complex<float> >(rows,cols,100, 0.001,"1.3","complex<float>","./cgfltcplx.txt",cout);
    CGTest< complex<double> >(rows,cols,100, 0.001,"1.4","complex<double>","./cgdblcplx.txt",cout);

  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
