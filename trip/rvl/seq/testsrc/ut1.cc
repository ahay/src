#include "seqfos.hh"

using namespace RVL;

int main() {

  SeqDC x;
  
  std::list<double> a;
  for (int i=0;i<6;i++) {
    a.push_back(double(i*i));
  }

  AssignList f(a);
  
  std::vector<DataContainer const *> v;
  x.eval(f,v);

  x.write(cout);

}


    
