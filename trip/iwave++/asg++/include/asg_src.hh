#ifndef __ASG_SRC
#define __ASG_SRC

#include "asg_headers.h"

#include "alg.hh"
#include "state.hh"
#include "seamx_headers.hh"

namespace ASG {

  using RVLAlg::Algorithm;
  using RVL::RVLException;
  using TSOpt::IWaveState;

  class ASGSource: public Algorithm {
    
  private:

    tracegeom & tg;

    POINTSRC * src;
    SAMPLER * arr;

    mutable RPNT scoord;
    mutable int istart;

    IWaveState & state;

    mutable bool ready;
    mutable bool ans;

    ASGSource();
    ASGSource(ASGSource const &);

  public: 
    
    ASGSource(IWaveState & _state, tracegeom & _tg);
    ~ASGSource();
    void init();
    void run();
    int getStartIndex();
    void fprint();

  };
}

#endif
