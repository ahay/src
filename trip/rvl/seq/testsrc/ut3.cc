#include "seqspace.hh"

using namespace RVL;

int main() {

  SeqSpace sp;
  
  std::list<double> a;
  for (int i=0;i<6;i++) {
    a.push_back(double(i*i));
  }

  AssignList f(a);
  
  Vector<double> v(sp);
  v.eval(f);

  v.write(cout);

}


    
