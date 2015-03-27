#ifndef __IW_SAMP__
#define __IW_SAMP__

#include "except.hh"
#include "parser.h"
#include "parserdecl.hh"
#include "iwave.h"
#include "grid.h"
#include "traceio.h"

namespace TSOpt {
  using RVL::parse;
  using RVL::valparse;
  using RVL::RVLException;
  using RVL::ProtectedDivision;

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
    int taperwidth;   // width of taper for each side for one shot in (# of traces)
    int timewidth;    // width of time to taper at the end of simulation in (ms)

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
	  
}

#endif
