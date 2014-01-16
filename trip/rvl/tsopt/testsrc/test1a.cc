#include "rns.hh"
#include "dtterm.hh"
#include "tests.hh"

int main(int argc, char ** argv) {

  using namespace TSOpt;
  using namespace RVL;
  using namespace RVLAlg;

  size_t n=1;
  float dt=0.1;
  bool verbose=true;
  
  // define Sim object, initialize control
  ContentPackage<float,size_t> c;
  c.initialize(n);
  for (int i=0;i<n;i++) c.getData()[i]=0.0;
  Dyn1 stp(c, dt, verbose);

  cout << "current time   = "; stp.getTime().write(cout);
  cout << "next time      = "; stp.getNextTime().write(cout);

}
    
    
