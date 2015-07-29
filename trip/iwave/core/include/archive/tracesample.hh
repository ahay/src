#ifndef __IWAVE_TRACE
#define __IWAVE_TRACE

#include "seamx_headers.hh"
#include "state.hh"
#include "alg.hh"

namespace TSOpt {

  using RVL::RVLException;
  using RVLAlg::Algorithm;

  /** object wrapper for trace term */
  template<typename State>
  class TraceSampler: public Algorithm {

  private:

    State & state;
    SAMPLER & trace;

    TraceSampler();
    TraceSampler(TraceSampler<State> const &);

  public:
    
    TraceSampler(State &_state, SAMPLER & _trace) 
      : state(_state), trace(_trace) {}
    virtual ~TraceSampler() {}
    //    SAMPLER & get_traceterm() {return trace;}
    //    SAMPLER const & get_traceterm() const {return trace;}
    virtual void run() {
      int err=0;
      if ((err=sampler_run(&trace,&(state.State::getIWAVE().model)))) {
	RVLException e;
	e<<"Error: TraceSampler::init from sampler_run\n";
	e<<"returned code "<<err<<"\n";
	throw e;
      }
    }
  
    ostream & write(ostream & str) const {
      str<<"TraceSampler: output to C stream\n";
      sampler_fprint(&trace,state.getStream());
      return str;
    }
    
  };

  /** object wrapper for trace term */
  template<typename State> 
  class NoTraceSampler: public TraceSampler<State> {

  private:

    NoTraceSampler();
    NoTraceSampler(NoTraceSampler<State> const &);

  public:
    
    NoTraceSampler(State &_state, SAMPLER & _trace)
      : TraceSampler<State>(_state,_trace) {}
    ~NoTraceSampler() {}
    void run() {}
    
    ostream & write(ostream & str) const {
      str<<"NoTraceSampler: vacuous trace sampler\n";
      return str;
    }

  };

}

#endif
