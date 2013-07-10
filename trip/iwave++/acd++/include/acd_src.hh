#ifndef __ASG_SRC
#define __ASG_SRC

#include "alg.hh"
#include "state.hh"
#include "seamx_headers.hh"
#include "acd_gfdm.h"

namespace ACD {

  using RVLAlg::Algorithm;
  using RVL::RVLException;
  using TSOpt::IWaveState;

  class ACDSource: public Algorithm {
    
  private:

    tracegeom & tg;
    SAMPLER * arr;
    mutable int istart;
    IWaveState & state;

    mutable bool ready;
    mutable bool ans;

    ACDSource();
    ACDSource(ACDSource const &);

  public: 
    
    ACDSource(IWaveState & _state, tracegeom & _tg);
    ~ACDSource();
    void init();
    void run();
    int getStartIndex();
    void fprint();

  };
}

#endif
