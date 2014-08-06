#include "rns.hh"
#include "dtterm.hh"
#include "tests.hh"
#include "sim.hh"
#include "ACPsim.hh"
#include "jet.hh"
#include "Dyn1_soln.hh"
#include "stdVector.hh"

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

  FwdDTimeTerm<RnState> tt(r);
  tt.setTargetTime(nt);
  
  initFwd<RnState> aa(stp);

  // make a sim
  ACPSim<RnState, StdVector<RnState, std::allocator<RnState> > > f(stp,tt,aa, 5);
  
  // stuff for linearized step
  ContentPackage<float,size_t> dc;
  dc.initialize(n);
  for (int i=0;i<n;i++) dc.getData()[i]=1.0;

  
  // Comment: dc is the some data container holding init. lin state
  //          r is the ref simulator's data container
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
  BwdDTimeTerm<RnState> att(ar);
  att.setTargetTime(0);

  initAdj<RnState> cc(astp);

  StdSim1Jet<RnState, RnState> j(f, dstp, bb, dtt, astp, cc, att);

  j.getAdj().run();

  // Check output

  cout << "Test Dyn1_Adj_Soln = " << Dyn1_Adj_Soln(0) << endl;
  cout << "astp.getState(0) = " <<  astp.getState().getrn().u[0] << endl; 
  
  
  if ( (Dyn1_Adj_Soln(0) - astp.getState().getrn().u[0]) < 1e-6) {
    return 0;
  }

  else {
    return 1;
  }
  

  return 0;

}
    
    
