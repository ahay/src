#include "rnsfo.hh"

namespace TSOpt {

  void LDCtoRnFOR::operator () (LocalDataContainer<float> const & x) {
    try {
      // initialize dimension - it=0
      // <ME?> what if s already exists?
      if (!( s.getrn().u && s.getrn().c && s.getrn().nu>=1 && s.getrn().nc>=1))  {
	s.initialize(x.getSize(),x.getSize());
      }
      
      memcpy(s.getrn().c,x.getData(),x.getSize()*sizeof(float));
      memcpy(s.getrn().u,x.getData(),x.getSize()*sizeof(float));

      StdDiscreteTime t;
      // <ME?> Is this incorrect? Should be at the initial time of 
      //       the simulation, whether fwd or adj
      t=0;
      s.setTime(t);
      
    
    }
    catch (RVLException & e) {
      e<<"\ncalled from LDCtoRnFOR::operator()\n";
      throw e;
    } 
  }

  void RntoLDCFO::operator()(LocalDataContainer<float> & x) {
    try {
      // check dimension
      if (s.getrn().nu != x.getSize()) {
	RVLException e;
	e<<"Error: RntoLDCFO::operator()\n";
	e<<"external, interal output dims don't agree\n";
	throw e;
      }
      out=x.getData();
    }
    catch (RVLException & e) {
      e<<"\ncalled from RntoLDCFO::operator()\n";
      throw e;
    } 
  }

  RntoLDCFO::~RntoLDCFO() {
    if (out) 
      memcpy(out,s.getrn().u,s.getrn().nu*sizeof(float));
  }

  /*
  void LDCtoRealRnFOR::operator () (LocalDataContainer<double> const & x) {
    try {
      // initialize dimension - it=0
      // <ME?> what if s already exists?
      if (!( s.getrealrn().u && s.getrealrn().c && s.getrealrn().nu>=1 && s.getrealrn().nc>=1))  {
	s.initialize(x.getSize(),x.getSize());
      }

      cout << "x.getData()[0]" << x.getData[0] << endl;
      

      int nc = s.getrealrn().nc;

      memcpy(s.getrealrn().c, c, nc*sizeof(double));
      memcpy(s.getrealrn().u, x.getData(), x.getSize()*sizeof(double));

      RealTime t;
      // <ME?> Is this incorrect? Should be at the initial time of 
      //       the simulation, whether fwd or adj
      t.assign(0.0,0);
      s.setTime(t);
      
    
    }
    catch (RVLException & e) {
      e<<"\ncalled from LDCtoRealRnFOR::operator()\n";
      throw e;
    } 
  }
  */

  void RealRntoLDCFO::operator()(LocalDataContainer<double> & x) {
    try {
      // check dimension
      if (s.getrealrn().nu != x.getSize()) {
	RVLException e;
	e<<"Error: RealRntoLDCFO::operator()\n";
	e<<"external, interal output dims don't agree\n";
	throw e;
      }
      out=x.getData();
    }
    catch (RVLException & e) {
      e<<"\ncalled from RealRntoLDCFO::operator()\n";
      throw e;
    } 
  }

  RealRntoLDCFO::~RealRntoLDCFO() {
    if (out) {
      memcpy(out,s.getrealrn().u,s.getrealrn().nu*sizeof(double));
    }
  }





}
