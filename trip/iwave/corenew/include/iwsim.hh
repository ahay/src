#ifndef __IW_STATE__
#define __IW_STATE__

#include "alg.hh"
#include "revolve.h"
#include "iwenv.hh"
#include "iwsamp.hh"
#include "iwtree.hh"

namespace TSOpt {
  using RVL::parse;
  using RVL::RVLException;
  using RVL::ProtectedDivision;
  using RVLAlg::Algorithm;

  class IWaveSim: public Algorithm {

  private:
    IWaveInfo const & ic;   // model defn struct
    IWaveTree * w;                   // wraps vector of IMODELs
    std::vector<TASK_RELN *> t;      // i/o tasks
    std::vector<IWaveSampler *> s;   // i/o samplers
    grid g;                          // simulation grid
    bool fwd;
    FILE * stream;

    // related to checkpoint alg
    int printact;
    int order;
    int snaps;
    int ndyn;
    RARR *** cps;
    int narr;

    bool dryrun;
    ostream & drystr;

    ostream & announce;

    IWaveSim();
    IWaveSim(IWaveSim const &);

  public:
    IWaveSim(int order, 
	     bool fwd, 
	     PARARRAY & par, 
	     FILE * stream, 
	     IWaveInfo const & _ic, 
	     int printact=0,
	     int snaps=0, 
	     bool dryrun=false,
	     ostream & drystr=cerr,
	     ostream & announce=cerr);
    ~IWaveSim();
    void run();
    std::vector<IWAVE *> const & getStateArray() const { return w->getStateArray(); }
    std::vector<RDOM *> const & getRDOMArray() const { return w->getRDOMArray(); }
    void printgrid(FILE * fp) const;
    ostream & write(ostream & str) const;
  };

  void IWaveApply(int argc, char ** argv);

}

#endif
