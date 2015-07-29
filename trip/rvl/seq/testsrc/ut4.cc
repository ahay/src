#include "multop.hh"

using namespace RVL;

int main() {

  SeqSpace sp;

  Vector<double> va(sp);
  Vector<double> vb(sp);
  Vector<double> vo(sp);
  
  std::list<double> a;
  for (int i=0;i<3;i++) {
    a.push_back(double(i));
  }

  std::list<double> b;
  for (int i=0;i<2;i++) {
    b.push_back(1.0);
  }

  AssignList fa(a);
  va.eval(fa);

  PolyMultOp op(b);

  op.applyOp(va,vo);

  cout<<"*********************************************\n";
  cout<<"*   Polynomial Multiplcation Operator Test  *\n";
  cout<<"*********************************************\n";
  cout<<"\ninput vector: \n";
  va.write(cout);
  cout<<"\nmultiplication operator \n";
  op.write(cout);
  cout<<"\noutput vector:\n";
  vo.write(cout);

}


    
