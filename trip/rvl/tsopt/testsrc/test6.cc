#include "rnspace.hh"
#include "rns.hh"
#include "dtterm.hh"
#include "tests.hh"
#include "rnsfo.hh"
#include "sim.hh"
#include "jet.hh"
#include "simop.hh"
#include "pol.hh"

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
  Dyn1 stp(dt, verbose);
  RnState & r = dynamic_cast<RnState &>(stp.getState());

  FwdDTimeTerm<RnState> tt(r);
  tt.setTargetTime(nt);
  
  initFwd_RVL<RnState> aa(stp);
  
  // make a sim
  StdSim<RnState> f(stp,tt,aa);

  // make 1-jet
  D_Dyn1 dstp(r, dt, verbose);
  RnState & dr = dynamic_cast<RnState &>(dstp.getState());
  FwdDTimeTerm<RnState> dtt(dr);
  dtt.setTargetTime(nt);
  initLin_RVL<RnState> bb(dstp);

  A_Dyn1 astp(r, dt, verbose); 
  RnState & ar = dynamic_cast<RnState &>(astp.getState());
  BwdDTimeTerm<RnState> att (ar);
  att.setTargetTime(0);
  initAdj_RVL<RnState> cc(astp);

  StdSim1Jet<RnState, RnState> j(f, dstp, bb, dtt, astp, cc,att);
  
  // make spaces
  RnSpace<float> dom(n);
  RnSpace<float> rng(n);

  // make vector input
  Vector<float> x(dom);
  x.zero();

  // make simop
  SimOp
    <
    float,  
    RnState,
    OpNewCreatePolicy<RnState,LDCtoRnFOR>,
    OpNewConstCreatePolicy<RnState,RntoLDCFO> 
    > op(dom,rng,j);

  OperatorEvaluation<float> opeval(op,x); 

  cout<<"***** TSOPT Test 6 *****"<<endl;
  cout<<"target result: length 1 vector with sample =  0.0897965"<<endl;
  cout<<"computed result:\n";
  opeval.getValue().write(cout);
  
}
    
    
