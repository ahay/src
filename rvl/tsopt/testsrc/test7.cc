

#include "rnspace.hh"
#include "rns.hh"
#include "dtterm.hh"
#include "tests.hh"
#include "rnsfo.hh"
#include "sim.hh"
#include "jet.hh"
#include "simop.hh"
#include "pol.hh"
#include "adjtest.hh"
#include "CPsim.hh"
#include "ACPsim.hh"
#include "stdVector.hh"

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

  // define forward step/term
  Dyn1 stp(dt, verbose);
  RnState & r = dynamic_cast<RnState &>(stp.getState());
  FwdDTimeTerm<RnState> tt(r);
  tt.setTargetTime(nt);

 
  NullStateAlg<RnState> NullInit(r);
  initFwd_RVL<RnState> aa(stp);
    
  // make a sim

  StdSim<RnState> f(stp,tt,aa);
  //CPSim<RnState, StdVector<RnState, std::allocator<RnState> > > f(stp,tt, aa, nt, 3);
  //RASim<RnState, StdVector<RnState, std::allocator<RnState> > > f(stp,tt, aa);
  //ACPSim<RnState, StdVector<RnState, std::allocator<RnState> > > f(stp,tt, aa, 5);
  
  // make derivative step / term
  D_Dyn1 dstp(r, dt, verbose);
  RnState & dr = dynamic_cast<RnState &>(dstp.getState());
  FwdDTimeTerm<RnState> dtt(dr);
  dtt.setTargetTime(nt);
  initLin_RVL<RnState> bb(dstp);

  // make adjoint step / term
  A_Dyn1 astp(r, dt, verbose); 
  RnState & ar = dynamic_cast<RnState &>(astp.getState());
  BwdDTimeTerm<RnState> att (ar);
  att.setTargetTime(0);
  initAdj_RVL<RnState> cc(astp);

  // make jet object
  StdSim1Jet<RnState, RnState> j(f, dstp, NullInit,  dtt, astp, NullInit,  att);
  
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

  RVLRandomize<float> adc(15061949,-1.0,1.0);
  
  // AdjointTest
  AdjointTest(opeval.getDeriv(),adc, cout);
 

 
}
    
    
