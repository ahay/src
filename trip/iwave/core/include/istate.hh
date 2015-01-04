#ifndef __IW_STATE__
#define __IW_STATE__

#include "alg.hh"
#include "except.hh"
#include "parser.h"
#include "parserdecl.hh"
#include "iwinfo.hh"
#include "iwave.h"
#include "grid.h"
#include "traceio.h"
#include "revolve.h"

namespace TSOpt {
  using RVL::parse;
  using RVL::valparse;
  using RVL::RVLException;
  using RVL::ProtectedDivision;
  using RVLAlg::Algorithm;

  /** Function to initialize parallel environment, output stream, and param table */
  void IWaveEnvironment(int argc, char ** argv, int ts, PARARRAY ** par, FILE ** stream);

  class IWaveSampler {
  private:
    string samplekey;
    string suffix;
    string pname;
    // FILE * stream; 
    std::vector<axis *> axes; 
    axis const * getTimeAxis(int dim) const;
    bool has_Axis(int i) const;
    bool has_Spatial_Axes;
    mutable int prev_panelindex;

    // this object really should be private data
    // of a trace sampler subclass or independent object
    tracegeom * tg;
    int sampord;
    int tracestart;
    int tracestop;
    int dump_term;

  public:
    IWaveSampler(IWAVE * state, string key, PARARRAY & pars, FILE * stream);
    virtual ~IWaveSampler();
    int getNumAxes() const { return axes.size(); }
    axis const & getAxis(int i) const;
    ireal getCellVol() const;
    ireal getRecipCellVol() const;
    void sample(grid g, IPNT step, bool fwd, bool input, 
		IWAVE * state, int ridx, int iwdx, FILE * stream,
		bool dryrun=false, ostream & drystr=cerr);
  };
	  
  // array of IWaveStates with ref object of same
  // type pointing to first half of array
  class IWaveTree {
  private:
    bool own;
    std::vector<IWAVE *> sa;    std::vector<RDOM *> rd;
    IWaveTree * ref;
    // required for proper destruction
    IWaveInfo const & ic;
    IWaveTree();
    IWaveTree(IWaveTree const &);
    IWaveTree(std::vector<IWAVE *> sv, IWaveInfo const & _ic);

  public:
    IWaveTree(PARARRAY & _pars, FILE * _stream, IWaveInfo const & _ic,
	      int order=0);
    ~IWaveTree();
    
    std::vector<IWAVE *> & getStateArray() { return sa; }
    std::vector<IWAVE *> const & getStateArray() const { return sa; }
    std::vector<IWAVE *> & getRefStateArray() { 
      if (ref) return ref->getStateArray();
      else {
	RVLException e;
	e<<"ERROR: IWaveTree::getRefStateArray()\n";
	e<<"  ref state not initialized - probably order = 0\n";
	e<<"  so no derivative so no reference state\n";
	throw e;
      }
    }
    std::vector<IWAVE *> const & getRefStateArray() const { 
      if (ref) {return ref->getStateArray(); }
      else {
	RVLException e;
	e<<"ERROR: IWaveTree::getRefStateArray()\n";
	e<<"  ref state not initialized - probably order = 0\n";
	e<<"  so no derivative so no reference state\n";
	throw e;
      }
    }
    std::vector<RDOM *> const & getRDOMArray() const { return rd; }
    std::vector<RDOM *> const & getRefRDOMArray() const { return ref->getRDOMArray(); }

    ostream & write(ostream & str) const {
      str<<"IWaveTree, length "<<sa.size()<<"\n";
      return str;
    }
  };

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
