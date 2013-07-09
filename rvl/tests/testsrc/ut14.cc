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
#include "op.hh"
#include "derivtest.hh"

using namespace RVL;

class f: public BinaryLocalFunctionObject<float> {
public:
  f() {}
  f(const f & ) {}
  ~f() {}
  using RVL::BinaryLocalEvaluation<float>::operator();
  void operator()(LocalDataContainer<float> & x,
		  LocalDataContainer<float> const & y) {
    int n=x.getSize();
    for (int i=0;i<n;i++) {
      x.getData()[i] = (y.getData()[i])*(y.getData()[i]*y.getData()[i]);
    }
  }
  string getName() const { return "f"; }
};

class df: public TernaryLocalFunctionObject<float> {
public:
  df() {}
  df(const f & ) {}
  ~df() {}
  using RVL::TernaryLocalEvaluation<float>::operator();
  void operator()(LocalDataContainer<float> & x,
		  LocalDataContainer<float> const & y,
		  LocalDataContainer<float> const & dy) {
    int n=x.getSize();
    for (int i=0;i<n;i++) {
      x.getData()[i] = 3.0*(y.getData()[i])*(y.getData()[i])*(dy.getData()[i]);
    }
  }
  string getName() const { return "df"; }
};

int main() {
  try {
    PlantSeeds(getpid());

    cout<<"/*************************************************"<<endl;
    cout<<" *            BEGIN RVL UNIT TEST 14             *"<<endl;
    cout<<" * Test of Operator - domain/range sanity test   *"<<endl;
    cout<<" *************************************************/"<<endl; 
    cout<<endl<<"1. construct 10D RnSpace, several FOs"<<endl<<endl;
    int n=10;
    int m=11;
    RnSpace<float> sp(n);
    f fff;
    df dfff;
    
    cout<<endl<<"2. construct OpFO"<<endl<<endl;
    OpFO<float> doit(sp,sp,fff,dfff,dfff);

    cout<<endl<<"3. construct different space, vector in it"<<endl<<endl;
    RnSpace<float> wsp(m);
    Vector<float> x(wsp);
    float c = 2.0;
    RVLAssignConst<float> ac(c);
    x.eval(ac);

    cout<<endl<<"4. attempt to evaluate - should throw exception"<<endl<<endl;
    OperatorEvaluation<float> doiteval(doit,x);

    cout<<endl<<"5. if you get this far, something is wrong"<<endl<<endl;

    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *              END RVL UNIT TEST 14             *"<<endl;
    cout<<" *************************************************/"<<endl; 

  }
  catch (RVLException & e) {
    ofstream str("testsrc/ut14/ut14.out");
    cout<<endl;
    cout<<"successful conclusion - error trapped - see output in ut14/ut14.out\n";
    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *              END RVL UNIT TEST 14             *"<<endl;
    cout<<" *************************************************/"<<endl;     
    e.write(str);
    exit(0);
  }
}
