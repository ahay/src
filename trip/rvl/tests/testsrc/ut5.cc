/*************************************************************************

Copyright Rice University, 2004.
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

#include "rnspace.hh"
#include "functional.hh"
#include "gradtest.hh"

using namespace RVL;

/** Test functional - computes square of square norm */

class mnsq: public Functional<float> {

private:

  RnSpace<float> & sp;
  mnsq();

protected:

  virtual Functional<float> * clone() const {
    return new mnsq(*this);
  }

  virtual Space<float> & getDomain() const { return sp; }

  virtual void apply(const Vector<float> & x, float & val) const {
    val = (x.inner(x))*(x.inner(x));
  }

  virtual void applyGradient(const Vector<float> & x, 
			Vector<float> & g) const {
    g.copy(x);
    g.scale(4.0*(x.inner(x)));
  }

  virtual void applyHessian(const Vector<float> & x,
			    const Vector<float> & yin, 
			    Vector<float> & yout) const {}

public:

  mnsq(const mnsq & fcn): sp(fcn.sp) {}
  mnsq(RnSpace<float> & _sp): sp(_sp) {}
  ~mnsq(){}

  ostream & write(ostream & str) const {
    str<<"mnsq"<<endl;
    return str;
  }
};

int main() {
  try {
    //    srand(getpid());
    //srand(19490615);
    PlantSeeds(19490615);
    int n=5;
    RVLRandomize<float> f;
    RnSpace<float> sp(n);
    Vector<float> v(sp);
    v.eval(f);
    mnsq op(sp);

    FunctionalEvaluation<float> feval(op,v);
    cout<<"/*************************************************"<<endl;
    cout<<" *            BEGIN RVL UNIT TEST 5              *"<<endl;
    cout<<" * Test of Functional, CheckGradient and other   *"<<endl;
    cout<<" * methods - Functional implements               *"<<endl;
    cout<<" *          f(x) = (|x|^2)^2                     *"<<endl;
    cout<<" *************************************************/"<<endl; 
    cout<<endl<<"1. construct 5D RnSpace, random input Vector x, write out:"<<endl<<endl;
    v.write(cout);
    cout<<endl<<"2. compute value of f(x) = "<<feval.getValue()<<endl;

    cout<<"   compare to result of independent computation\n";
    cout<<"   using RVL::Vector::inner = "<<(v.inner(v))*(v.inner(v))<<endl;
    Vector<float> p(sp);
    p.eval(f);
    cout<<endl<<"3. construct random direction vector:"<<endl;
    p.write(cout);
    cout<<endl<<"4. check gradient against finite difference approximation\n";
    cout<<"   of directional derivative in direction p, using checkGradient.\n";
    cout<<"   output: last column should converge to 2.0, indicating 2nd\n";
    cout<<"   order convergence of centered diff. approx. to dir.  deriv.,\n";
    cout<<"   computed as inner product of gradient and direction."<<endl;
    cout<<"   convergence will start to degrade as relative error gets close\n";
    cout<<"   to level of 32 bit machine precision, i.e. 10^{-5 or 6}."<<endl;
    cout<<"   precise specs: success if convergence rate estimate >= 1.95"<<endl;
    cout<<"   test history reported in ut5.aux"<<endl<<endl;
    ofstream str("testsrc/ut5/ut5.out");
    if (!(GradientTest<float>(op,v,p,str,11,0.01,1.0))) {
      RVLException e;
      e<<"Error: ut5\n";
      e<<"gradient check failed\n";
      throw e;
    }
    else {
      cout<<"GradientTest returned true - success"<<endl;
    }
    cout<<"--- see report in testsrc/ut5/ut5.out"<<endl;
    str.close();
    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *              END RVL UNIT TEST 5              *"<<endl;
    cout<<" *************************************************/"<<endl; 
  }
  catch (RVLException & e) {
    e.write(cout);
    exit(1);
  }
}

  







