#include "seqfos.hh"

using namespace RVL;

int main() {

  SeqDC pa;
  SeqDC pb;
  SeqDC q;
  
  std::list<double> a;
  for (int i=0;i<3;i++) {
    a.push_back(double(i));
  }

  std::list<double> b;
  for (int i=0;i<2;i++) {
    b.push_back(1.0);
  }

  AssignList fa(a);
  AssignList fb(b);

  std::vector<DataContainer const *> v;
  pa.eval(fa,v);
  pb.eval(fb,v);

  PolyMult m(pb);
  v.push_back(&pa);
  q.eval(m,v);

  cout<<"*********************************************\n";
  cout<<"*       Polynomial Multiplcation Test       *\n";
  cout<<"*********************************************\n";
  cout<<"\nfirst factor: \n";
  pa.write(cout);
  cout<<"\nsecond factor \n";
  pb.write(cout);
  cout<<"\nproduct:\n";
  q.write(cout);

}


    
