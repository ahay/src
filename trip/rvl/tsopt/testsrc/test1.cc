#include "cdt.hh"

int main(int argc, char ** argv) {

  using namespace TSOpt;
  using namespace RVL;

  cout<<"create DiscreteTime object\n";  
  StdDiscreteTime t1;
  cout<<"assign zero, check\n";
  t1=0;
  t1.write(cout);
  cout<<"create another discrete time object\n";
  StdDiscreteTime t2;
  cout<<"assign it a different integer\n";
  t2=11;
  cout<<"assign second DT to first\n";
  t1=t2;
  cout<<"write them both out\n";
  t1.write(cout);
  t2.write(cout);

}
    
    
