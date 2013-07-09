

#include "rnspace.hh"
#include "rns.hh"
#include "dtterm.hh"
#include "tests.hh"
#include "rnsfo.hh"
#include "sim.hh"
#include "jet.hh"
#include "simop.hh"
#include "pol.hh"
#include "derivtest.hh"


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
  
  Dyn1 stp(dt, verbose);
  RnState & r = dynamic_cast<RnState &>(stp.getState());

  FwdDTimeTerm<RnState> tt(r);
  tt.setTargetTime(nt);

  NullStateAlg<RnState> NullInit(r);
  initFwd_RVL<RnState> aa(stp);
  
  // make a sim
  StdSim<RnState> f(stp,tt,NullInit);

  // make 1-jet
  D_Dyn1 dstp(r, dt, verbose);
  RnState & dr = dynamic_cast<RnState &>(dstp.getState());
  FwdDTimeTerm<RnState> dtt(dr);
  dtt.setTargetTime(nt);
  initLin_RVL<RnState> bb(dstp);

  verbose = true;
  A_Dyn1 astp(r, dt, verbose); 
  RnState & ar = dynamic_cast<RnState &>(astp.getState());
  BwdDTimeTerm<RnState> att (ar);
  att.setTargetTime(0);
  initAdj_RVL<RnState> cc(astp);
 
  StdSim1Jet<RnState, RnState> j(f, dstp, NullInit, dtt, astp, NullInit, att);
  
  // make spaces
  RnSpace<float> dom(n);
  RnSpace<float> rng(n);

 
  // make simop
  SimOp
    <
    float,  
    RnState,
    OpNewCreatePolicy<RnState,LDCtoRnFOR>,
    OpNewConstCreatePolicy<RnState,RntoLDCFO> 
    > op(dom,rng,j);


  // make vector inputs
  Vector<float> c(dom);
  Vector<float> dc(dom);

  RVLRandomize<float> adc(15061949,-1.0,1.0);
  dc.eval(adc);
  
#ifdef MAXPATHLEN
#define PATHLEN MAXPATHLEN
#else
#define PATHLEN 1024
#endif
  char buf[PATHLEN];
  memset(buf,0,PATHLEN);
  getcwd(buf,PATHLEN);
  std::string pth = buf;
  std::string pn = pth + "/derivtest.rpt";
  ofstream str(pn.c_str());
  cout<<"DerivTest returns "<<DerivTest(op, c, dc, str)<<"; full output in derivtest.rpt";





}
    
    
