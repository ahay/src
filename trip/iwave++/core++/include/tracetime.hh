#ifndef __IWAVE_TRACETIME
#define __IWAVE_TRACETIME

#include "state.hh"
#include "tterm.hh"

namespace TSOpt {

  template<typename State=IWaveState>
  class FwdTraceTime: public FwdTimeTerm<State> {

  private:

    TSIndex tgt;
    TIMESTEPINDEX _t;
    FwdTraceTime();
    FwdTraceTime(FwdTraceTime<State> const &);

  public:
    
    FwdTraceTime(State const & state,
		 bool verbose=false,
		 ostream & str=cout)
      : FwdTimeTerm<State>(state,verbose,str), tgt(_t) {}

    void setTargetTime(Time const & _t0) {
      try {
	tgt=_t0;
      }
      catch (bad_cast) {
	RVLException e;
	e<<"asdioiuthpi\n";
	throw e;
      }
    }
    Time const & getTargetTime() const { return tgt; }

    /** for use during initialization */
    void setTargetTime(int nt, float dt) { 
      _t.it=nt; _t.iv=0; _t.dt=dt; 
    }
    
  };

  template<typename State> 
  class BwdTraceTime: public BwdTimeTerm<State> {

  private:

    TIMESTEPINDEX _t;
    TSIndex tgt;
    BwdTraceTime();
    BwdTraceTime(BwdTraceTime<State> const &);

  public:
    
    BwdTraceTime(State const & state,
		 bool verbose=false,
		 ostream & str=cout)
      : BwdTimeTerm<State>(state,verbose,str), tgt(_t) {}

    void setTargetTime(Time const & _t0) { tgt=_t0; }
    Time const & getTargetTime() const { return tgt; }

    /** for use during initialization */
    void setTargetTime(int nt, float dt) { _t.it=nt; _t.iv=0; _t.dt=dt; }

    /** new query function (it <= itstop) **/
    //     bool query() { 
    //       try {
    // 	if (this->getTargetTime() >  this->s.State::getTime()) {
    // 	  return true; 
    // 	}
    // 	return false;
    //       }
    //       catch (RVLException & e) {
    // 	e<<"\n called from BwdTraceTerm::query\n";
    // 	throw e;
    //       }
    //     }
    
  };

}

#endif
