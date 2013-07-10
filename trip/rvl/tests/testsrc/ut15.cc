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

class a: public BinaryLocalFunctionObject<float> {
public:
  a() {}
  a(const a &) {}
  ~a() {}
  using RVL::BinaryLocalEvaluation<float>::operator();
  void operator()(LocalDataContainer<float> & x,
		  LocalDataContainer<float> const & y) {
    int n=x.getSize();
    for (int i=0;i<n;i++) {
      x.getData()[i] = float(i+1)*y.getData()[i];
    }
  }
  string getName() const { return "a"; }
};

class g: public BinaryLocalFunctionObject<float> {
public:
  using RVL::BinaryLocalEvaluation<float>::operator();
  void operator()(LocalDataContainer<float> & x,
		  LocalDataContainer<float> const & y) {
    int n=x.getSize();
    int m=y.getSize();
    if (n != 2) {
      RVLException e; 
      e<<"boo boo g\n";
      throw e;
    }
    if (m != 3) {
      RVLException e; 
      e<<"boo boo g\n";
      throw e;
    }
    for (int i=0;i<2;i++) {
      x.getData()[i] = 1.0f + y.getData()[i]-y.getData()[2];
    }
  }
  string getName() const { return "g"; }
};

class dg: public TernaryLocalFunctionObject<float> {
public:
  using RVL::TernaryLocalEvaluation<float>::operator();
  void operator()(LocalDataContainer<float> & x,
		  LocalDataContainer<float> const & y,
		  LocalDataContainer<float> const & dy) {
    int n=x.getSize();
    int m=y.getSize();
    if (n != 2) {
      RVLException e; 
      e<<"boo boo g\n";
      throw e;
    }
    if (m != 3) {
      RVLException e; 
      e<<"boo boo g\n";
      throw e;
    }
    for (int i=0;i<2;i++) {
      x.getData()[i] = dy.getData()[i]-dy.getData()[2];
    }
  }
  string getName() const { return "dg"; }
};

class ag: public TernaryLocalFunctionObject<float> {
public:
  using RVL::TernaryLocalEvaluation<float>::operator();
  void operator()(LocalDataContainer<float> & dy,
		  LocalDataContainer<float> const & y,
		  LocalDataContainer<float> const & x) {
    int n=x.getSize();
    int m=y.getSize();
    if (n != 2) {
      RVLException e; 
      e<<"boo boo ag\n";
      throw e;
    }
    if (m != 3) {
      RVLException e; 
      e<<"boo boo ag\n";
      throw e;
    }
    dy.getData()[2]=0;
    for (int i=0;i<2;i++) {
      dy.getData()[2] += x.getData()[i];
      dy.getData()[i] = x.getData()[i];
    }
  }
  string getName() const { return "ag"; }
};
    
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
    cout<<" *            BEGIN RVL UNIT TEST 15             *"<<endl;
    cout<<" * Test of OpComp - compose several operators    *"<<endl;
    cout<<" *************************************************/"<<endl; 
    cout<<endl<<"1. construct 10D RnSpace, several FOs"<<endl<<endl;
    int n=3;
    int m=2;
    RnSpace<float> sp(n);
    RnSpace<float> sp1(m);
    a aaa;
    g ggg;
    dg dggg;
    ag aggg;
    f fff;
    df dfff;
    
    cout<<endl<<"2. construct OpFOs"<<endl<<endl;
    OpFO<float> doit(sp,sp,fff,dfff,dfff);
    OpFO<float> doitagain(sp,sp1,ggg,dggg,aggg);
    LinearOpFO<float> yetagain(sp1,sp1,aaa,aaa);

    cout<<endl<<"3. construct Op comp"<<endl<<endl;
    /*
    OpComp<float> oc;
    oc.setNext(doit);
    oc.setNext(doitagain);
    oc.setNext(yetagain);
    */
    OpComp<float> oc(doit, doitagain, yetagain);

    cout<<endl<<"4. invent input vector, temporaries\n";
    LocalVector<float> x(sp);
    x.getData()[0]=1.0f;
    x.getData()[1]=2.0f;
    x.getData()[2]=3.0f;
    x.write(cout);

    cout<<endl<<"5. evaluate"<<endl<<endl;
    OperatorEvaluation<float> oceval(oc,x);
    OperatorEvaluation<float> doiteval(doit,x);
    OperatorEvaluation<float> doitagaineval(doitagain,doiteval.getValue());
    OperatorEvaluation<float> yetagaineval(yetagain,doitagaineval.getValue());

    cout<<endl<<" evaluation by hand:\n";
    yetagaineval.getValue().write(cout);
    cout<<endl<<" evaluation by OpComp:\n";
    oceval.getValue().write(cout);

    cout<<endl<<"6. evaluate derivative, in direction x\n";
    cout<<endl<<" evaluation by hand:\n";
    LocalVector<float> x1(sp);
    LocalVector<float> y1(sp1);
    LocalVector<float> y(sp1);
    LocalVector<float> xa1(sp);
    LocalVector<float> ya1(sp1);
    LocalVector<float> xa(sp);
    doiteval.getDeriv().applyOp(x,x1);
    doitagaineval.getDeriv().applyOp(x1,y1);
    yetagaineval.getDeriv().applyOp(y1,y);
    y.write(cout);
    cout<<endl<<" evaluation by OpComp:\n";
    oceval.getDeriv().applyOp(x,y);
    y.write(cout);

    cout<<endl<<"7. evaluate adjoint derivative, in direction of deriv output\n";
    cout<<endl<<" evaluation by hand\n";
    yetagaineval.getDeriv().applyAdjOp(y,ya1);
    /*
    cout<<"hand adj step 2:\n";
    cout<<"input = \n";
    y.write(cout);
    cout<<"output = \n";
    ya1.write(cout);
    */
    doitagaineval.getDeriv().applyAdjOp(ya1,xa1);
    /*
    cout<<"hand adj step 1:\n";
    cout<<"input = \n";
    ya1.write(cout);
    cout<<"output = \n";
    xa1.write(cout);
    */
    doiteval.getDeriv().applyAdjOp(xa1,xa);
    /*
    cout<<"hand adj step 0:\n";
    cout<<"input = \n";
    xa1.write(cout);
    cout<<"output = \n";
    */
    xa.write(cout);
    cout<<endl<<" evaluation by OpComp:\n";
    oceval.getDeriv().applyAdjOp(y,xa);
    xa.write(cout);

    cout<<endl;
    cout<<"/*************************************************"<<endl;
    cout<<" *              END RVL UNIT TEST 15             *"<<endl;
    cout<<" *************************************************/"<<endl; 

  }
  catch (RVLException & e) {
    e.write(cerr);
    exit(1);
  }
}
