#include "rns.hh"
#include "dtterm.hh"
#include "tests.hh"
#include "sim.hh"
#include "jet.hh"
#include "Dyn1_soln.hh"

int main(int argc, char ** argv) {

  using namespace TSOpt;
  using namespace RVL;
  using namespace RVLAlg;

  int n=1;
  int nt=9;
  StdDiscreteTime ft;
  float dt=0.01;
  bool verbose=false;
  
  ft=nt;

  // define Sim object, initialize control
  ContentPackage<float,size_t> c;
  c.initialize(n);
  for (int i=0;i<n;i++) c.getData()[i]=0.5;
  Dyn1 stp(c, dt, verbose);
  RnState & r = dynamic_cast<RnState &>(stp.getState());

  
  initFwd<RnState> aa(stp);

  FwdDTimeTerm<RnState> tt(r);
  tt.setTargetTime(nt);
  
  // make a sim
  StdSim<RnState> f(stp,tt,aa);

  // make a 1-jet
  ContentPackage<float,size_t> dc;
  dc.initialize(n);
  for (int i=0;i<n;i++) dc.getData()[i]=1.0;
  D_Dyn1 dstp(dc, r, dt, verbose);
  RnState & dr = dynamic_cast<RnState &>(dstp.getState());
  FwdDTimeTerm<RnState> dtt(dr);
  dtt.setTargetTime(nt);

  initLin<RnState> bb(dstp);

  // stuff for adjoint step
  ContentPackage<float,size_t> ac;
  ac.initialize(n);
  for (int i=0;i<n;i++) ac.getData()[i]=0.5;

  A_Dyn1 astp(ac, r, dt, verbose); 
  RnState & ar = dynamic_cast<RnState &>(astp.getState());
  FwdDTimeTerm<RnState> att(ar);
  att.setTargetTime(0);

  initAdj<RnState> cc(astp);

  StdSim1Jet<RnState, RnState> j(f,dstp, bb ,dtt, astp, cc, att);

  j.getFwd().run();



  // Check output
  
   cout << "Test Dyn1_Fwd_Soln = " << Dyn1_Fwd_Soln(nt) << endl;
  cout << "stp.getState(9) = " <<  stp.getState().getrn().u[0] << endl; 

  /*
  float eps = 100*numeric_limits<float>::epsilon();
  if ( (Dyn1_Fwd_Soln(nt) - stp.getState().getrn().u[0]) < eps) {
    return 0;
  }
  

  else {
    return 1;
  }
  */

  return 0;
}
    
    
