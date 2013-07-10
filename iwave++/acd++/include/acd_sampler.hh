#ifndef __ACD_SAMP
#define __ACD_SAMP

#include "parserdecl.hh"
#include "coords.hh"
#include "samp.hh"
#include "tsind.hh"
#include "tracesample.hh"
#include "tracetime.hh"
#include "seamx_headers.hh"
#include "acd_src.hh"

namespace ACD {

  // from tsopt
  using TSOpt::Time;
  using TSOpt::Sampler;
  using TSOpt::LinSampler;
  using TSOpt::Coords;
  using TSOpt::TimeTerm;

  // from iwave++
  using TSOpt::TSIndex;
  using TSOpt::IWaveState;
  using TSOpt::IWaveLinState;
  using TSOpt::TraceSampler;
  using TSOpt::NoTraceSampler;
  using TSOpt::FwdTraceTime;
  using TSOpt::BwdTraceTime;
  
  // from alg
  using RVLAlg::Terminator;
  using RVLAlg::Algorithm;
  using RVLAlg::NoAlg;
  
  // from rvl
  using RVL::RVLException;
  using RVL::parse;

  class ACDSampler: public Sampler {
    
    IWaveState & state;
    mutable SAMPLER trace;
    TIMESTEPINDEX ct, cte;
    mutable FwdTraceTime<> tt;
    mutable TSIndex tstart;
    mutable TSIndex tfinal;
    mutable ACDSource aps;
    mutable TSOpt::TraceSampler<IWaveState> trs;
    mutable bool is_samplerinit; /**< indicates if samplerinit called for current record */

    ACDSampler();
    ACDSampler(ACDSampler const &);

  public:

    ACDSampler(IWaveState & _state);
    ~ACDSampler();
    /** inherited from sampler */
    bool init();
    /** to be done before time index updated, but AFTER state update */
    Algorithm & pre()  const { return aps; }
    /** to be done AFTER time index update AND AFTER state update */
    Algorithm & post() const { return trs; }
    TimeTerm & term() const { return tt; }
    Time const & start() const { return tstart; }
    Time const & end() const { return tfinal;}

    /** returns number of records in trace data */
    int getNumRecs() const;
    /** returns current record number */
    int getRecIndx() const;
    /** returns first record number in current group */
    int getFirstRec() const;
    /** returns last record number in current group */
    int getLastRec() const;

    /** flush output data to disk */
    void flush() const;
    /** since iwave output is to C stdio, and it's not
	easily possible to interchange with std::stream i/o,
	this method is a fake - does not actually write to 
	stream, writes instead to the state's FILE * */
    ostream & write(ostream & str) const;
  };

  class ACDLinSampler: public LinSampler {

  private:
    IWaveLinState & state;
    mutable SAMPLER trace;
    TIMESTEPINDEX ct, cte;
    mutable FwdTraceTime<> tt;
    mutable TSIndex tstart;
    mutable TSIndex tfinal;
    mutable FwdTraceTime<IWaveLinState> ltt;
    mutable ACDSource aps;
    mutable TSOpt::NoTraceSampler<IWaveState> trs;
    mutable TSOpt::TraceSampler<IWaveLinState> ltrs;
    mutable NoAlg abs; // formerly ACD_Bornsrc
    mutable bool is_samplerinit;     /**< indicates if samplerinit called for current record */    

    ACDLinSampler();
    ACDLinSampler(ACDLinSampler const &);
  public:
    ACDLinSampler(IWaveLinState & _state);
    ~ACDLinSampler();
    bool init();
    // NOTE: pre and post now mean PRE- and POST- time index update, 
    // ALL to be done after state updates
    Algorithm & pre()  const { return aps; }
    Algorithm & post() const { return trs; }
    Algorithm & linpre() const { return abs; } // this is now a no-op
    Algorithm & linpost() const { return ltrs; }
    TimeTerm & term() const { return tt; }
    TimeTerm & linterm() const {return ltt;}
    Time const & start() const { return tstart; }
    Time const & end() const { return tfinal;}

    /** returns number of records in trace data */
    int getNumRecs() const;
    /** returns current record number */
    int getRecIndx() const;
    /** returns first record number in current group */
    int getFirstRec() const;
    /** returns last record number in current group */
    int getLastRec() const;

    /** flush output data to disk */
    void flush() const;
    /** since iwave output is to C stdio, and it's not
	easily possible to interchange with std::stream i/o,
	this method is a fake - does not actually write to 
	stream, writes instead to the state's FILE * */
    ostream & write(ostream & str) const;
  };

  class ACDAdjSampler: public LinSampler {

    IWaveLinState & state;
    mutable SAMPLER trace;
    TIMESTEPINDEX ct, cte;

    mutable FwdTraceTime<> tt;
    mutable BwdTraceTime<IWaveLinState> bwdtt;
    mutable TSIndex tstart;
    mutable TSIndex tfinal;

    mutable ACDSource aps;
    /*dum-sampler (don't need output forward (ref) field) */
    mutable TSOpt::NoTraceSampler<IWaveState> trs;
    /* need adj action of TraceSampler to inject source for backward sim */
    mutable TSOpt::TraceSampler<IWaveLinState> adjtrs;
    mutable NoAlg abs; 
    mutable bool is_samplerinit;     /**< indicates if samplerinit called for current record */
    
    ACDAdjSampler();
    ACDAdjSampler(ACDAdjSampler const &);
  public:
    ACDAdjSampler(IWaveLinState & _state);
    ~ACDAdjSampler();
    bool init();
    Algorithm & pre()  const { return aps; }
    Algorithm & post() const { return trs; }
    /* imaging cond (adjoint born src)*/
    Algorithm & linpre() const { 
      return adjtrs;
      //return abs;
    }
    /* backword source (data residual) */
    Algorithm & linpost() const { 
       return abs; 
       //  return adjtrs; 
    }
    TimeTerm & term() const { return tt; }
    TimeTerm & linterm() const {return bwdtt;}
    Time const & start() const { return tstart; }
    Time const & end() const { return tfinal;}

    /** returns number of records in trace data */
    int getNumRecs() const;
    /** returns current record number */
    int getRecIndx() const;
    /** returns first record number in current group */
    int getFirstRec() const;
    /** returns last record number in current group */
    int getLastRec() const;

    /** flush output data to disk */
    void flush() const;

    /** since iwave output is to C stdio, and it's not
	easily possible to interchange with std::stream i/o,
	this method is a fake - does not actually write to 
	stream, writes instead to the state's FILE * */
    ostream & write(ostream & str) const;
  };

}
#endif


