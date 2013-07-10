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
#include "cgnealg.hh"
#include "LinFitLS.hh"

using namespace RVL;
using namespace RVLAlg;
using namespace RVLUmin;

template<typename T>
void LinFitTest(int rows, int cols, typename ScalarFieldTraits<T>::AbsType _nlvl, 
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
    str<<" * Purpose: test LinFitLS comp type ="<<datatype<<endl;
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
    Vector<T> err(rng);
    err.eval(rnd);
    T nlvl = _nlvl*rhs.norm()/err.norm();
    rhs.linComb(nlvl,err);
    str<<endl;

    str<<"4. Build, run CGNE algorithm \n";
    str<<endl;    
    strout<<endl<<"*******************************************************"<<endl;
    strout<<"* Test of LinFitLS Function"<<endl;
    strout<<"* max iterations       = "<<maxit<<endl;
    strout<<"* residual tolerance   = "<<rtol<<endl;
    strout<<"* normal res tolerance = "<<nrtol<<endl;
    strout<<"* trust radius         = "<<Delta<<endl;


    CGNEPolicyData<T> params;
    params.rtol=rtol;
    params.nrtol=nrtol;
    params.maxcount=maxit;
    params.verbose=false;

    // allocate rnorm, nrnorm scalars
    atype rnorm = numeric_limits<atype>::max();
    atype nrnorm = numeric_limits<atype>::max();

    // construct CGNEAlg
    //    CGNEAlg<T> iter(est,A,rhs,rnorm,nrnorm,rtol,nrtol,maxit,Delta,strout);
    CGNEPolicy<T> p;
    p.assign(params);
    CGNEAlg<T> * iter = p.build(est,A,rhs,rnorm,nrnorm,strout);
    
    // run it
    //    iter.run();
    iter->run();

      // Generate identity preconditioner
      GenMat<T> B(dom,dom);
      for (int i=0;i<cols;i++)
          B.getElement(i,i)=ScalarFieldTraits<T>::One();
    // construct LinFitLS function
    LinFitLS<T, CGNEPolicy<T>, CGNEPolicyData<T> > f(A,B,rhs,params,strout);
    est.zero();
    FunctionalEvaluation<T> feval(f,est);
    strout<<"value from LinFitLS    = "<<feval.getValue()<<endl;
    strout<<"rnorm from CGNE        = "<<rnorm<<endl;
    strout<<"rnorm from LinFitLS    = "<<sqrt(2.0*abs(feval.getValue()))<<endl;
    strout<<"difference rnorms      = "<<rnorm-sqrt(2.0*abs(feval.getValue()))<<endl;
    strout<<"*******************************************************"<<endl;
    strout.close();
    delete iter;
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

    LinFitTest<float>(rows,cols,0.0,100,0.001,0.001,1000000.0,"2.1","float","./cgneflt.rpt",cout);
 //   LinFitTest<double>(rows,cols,0.1,100,0.001,0.001,1000000.0,"2.2","double","./cgnedbl.rpt",cout);
//    LinFitTest<complex<float> >(rows,cols,0.1,100,0.001,0.001,1000000.0,"2.3","complex<float>","./cgnefltcplx.rpt",cout);
//    LinFitTest<double>(rows,cols,0.2,100,0.001,0.001,1000000.0,"2.4","complex<double>","./cgnedblcplx.rpt",cout);
//    LinFitTest<float>(rows,cols,0.0,100,0.001,0.001,1.0,"2.5","float","./cgneflt1.rpt",cout);
//    LinFitTest<double>(rows,cols,0.2,100,0.001,0.001,1.0,"2.6","complex<double>","./cgnedblcplx1.rpt",cout);
  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}

