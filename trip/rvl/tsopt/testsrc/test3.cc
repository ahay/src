#include "rns.hh"
#include "dtterm.hh"
#include "tests.hh"
#include "sim.hh"
#include "Dyn1_soln.hh"

int main(int argc, char ** argv) {

  using namespace TSOpt;
  using namespace RVL;
  using namespace RVLAlg;

  int n=1;
  int nt=9;
  float dt=0.01;
  bool verbose=false;
  
  // define Sim object, initialize control
  ContentPackage<float,size_t> c;
  c.initialize(n);
  for (int i=0;i<n;i++) c.getData()[i]=0.5;

  Dyn1 stp(c, dt, verbose);
  RnState & r = dynamic_cast<RnState &>(stp.getState());
  
   // final time terminator - set final time index
  FwdDTimeTerm<RnState> tt(r);
  tt.setTargetTime(nt);
  
 

  initFwd<RnState> aa(stp);
 
  // make a sim, run it
  StdSim<RnState> f(stp,tt,aa);
  f.run();


  // Check output

  cout << "Test Dyn1_Fwd_Soln = " << Dyn1_Fwd_Soln(nt) << endl;
  cout << "stp.getState(9) = " <<  stp.getState().getrn().u[0] << endl; 

  float eps = 100*numeric_limits<float>::epsilon();
  
  if ( (Dyn1_Fwd_Soln(nt) - stp.getState().getrn().u[0]) < eps) {
    return 0;
  }

  else {
    return 1;
  }

}
    
    
