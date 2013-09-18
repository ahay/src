#include "except.hh"
#include <deque>
#include "sim.hh"
#include "iwavetime.hh"

namespace TSOpt {

  using RVL::RVLException;

  size_t pow2(int);
  deque<bool> d2b(int);
  int b2d(deque<bool>);
  int level(deque<bool>);
  bool index_incr(deque<bool>, vector<int> &);
  void print_xcode(deque<bool> xcode, ostream & str);
  void print_idx(vector<int> idx, ostream & str);

  // array of IWaveStates with ref object of same
  // type pointing to first half of array
  class IWaveTreeState {
  private:
    bool own;
    std::vector<IWaveState *> sa;
    IWaveTreeState * ref;
    IWaveTreeState();
    IWaveTreeState(IWaveTreeState const &);
    IWaveTreeState(std::vector<IWaveState *> sv)
      : own(false), sa(iwave_max(sv.size()/2,0)) {
      if (sv.size()>1) 
	for (int i=0;i<sv.size()/2;i++) sa[i]=sv[i];
      ref=NULL;
    }
  public:
    IWaveTreeState(PARARRAY & _pars, FILE * _stream,
		   GFDM_INIT_FUN _minit, int n=0): own(true), sa(pow2(n)) {
      try {
	for (size_t i=0; i<sa.size(); i++) 
	  sa[i]=new IWaveState(_pars,_stream,_minit);
	// create reference state out of first half of
	// array, without allocating any new memory,
	// using private constructor
	if (n>0) ref = new IWaveTreeState(sa);
	else ref=NULL;
      }
    }
    ~IWaveTreeState() {
      if (ref) delete ref;
      if (own) 
	for (size_t i=0; i<sa.size(); i++) delete sa[i];
    }

    // use reference state (index 0) for all time references
    void setTime(Time const & t) {
      try { sa[0]->setTime(t);
      }
      catch (RVLException & e) {
	e<<"\ncalled from IWaveTreeState::setTime\n";
	throw e;
      }
    }
    Time const & getTime() const { return sa[0]->getTime(); }
    Time & getTime() { return sa[0]->getTime(); }
    Time & getNextTime() const { sa[0]->getTime(); }

    virtual std::vector<T *> & getStateArray() { return sa; }
    virtual std::vector<T *> const & getStateArray() const { return sa; }

    ostream & write(ostream & str) const {
      str<<"IWaveTreeState, length "<<sa.size()<<"\n";
      for (size_t i=0; i<sa.size(); i++) {
	str<<"** component "<<i<<":\n";
	sa[i]->write(str);
      }
      return str;
    }
  };

  // revise GFDM defn to include typedef void (*f)(std::vector<IWaveState> *> GFDSTEP;

  // Note: rename Mario's multi-sampler to trace2file - adjoint is file2trace
  IWaveFwdSim: public Algorithm {

  private:
    int order;                                    // derivative order
    IWaveTreeState iwstate;                       // wraps vector of IWaveStates
    GFD_MODEL gfdm;                               // forward time step function
    file2grid * rsf;                              // grid input
    file2trace * m_s_src;                         // trace input - for source
    trace2file * m_s_fwd;                         // trace output
    PARARRAY & pars;
    FILE * stream;

  public:
    IWaveFwdSim(int order, PARARRAY & _pars, FILE * _stream, GFD_INIT_FUN minit)
      : iwstate(_pars,_stream,minit,order), pars(_pars), stream(_stream) {
      minit(&gfdm);
      rsf = new file2grid;
      // need to think up constructor for this thing - combine possible shift
      // to dual grid with rsfread.
      // - use isdyn and specs->fd_set_grid_type to param file2grid
      m_s_src = new file2trace;
      if (file2trace_construct(m_s_src,&pars,hdrkey,stream)) {
	// throw exception
      }
      m_s_fwd = new trace2file;
      if (trace2file_construct(m_s_fwd,&pars,hdrkey,stream)) {
	// throw exception
      }
      // probably should be a movie option
    }

    void run() {
      try {

	while (ms_s_fwd->xrec <= ms_s_fwd->last) {
	  // zero everything
	  iwstate.zero();
	  // load sim data - used to be iwave_static_init
	  file2grid_run(rsf,iwstate,ms_s_fwd->xrec,par,stream);
	  // initialize sampler - always last component sampled
	  std::vector<IWaveState *> states = iwstate.getStateArray();
	  trace2file_init(ms_s_fwd,
			  &(states[states.size()-1]->model),
			  &(iwstate.getPARARRAY()),iwstate.getFILE());
	  // initialize source - always first component sampled
	  file2trace_init(ms_s_src,
			  &(states[0]->model),
			  &(iwstate.getPARARRAY()),iwstate.getFILE());
	    
	  // time loop - note that this requres storagee of the sampler
	  // start and stop times in TIMESTEPINDEX form
	  iwstate.setTime(ms_s_src.start);
	  while (less_than(iwstate.getTime(),ms_s_fwd.stop)) {
	    // time step
	    tsf(states);
	    // source input
	    file2trace_run(ms_s_src,&(states[0]->model));
	    // trace output
	    trace2file_run(ms_s_fwd,&(states[states.size()-1]->model));
	    // next time
	    next_step(iwstate.getTime());
	  } //end time loop 
	  
	  // write ouput - note that d, og inputs should be struct members
	  // and initialized by the _init function
	  trace2file_flush(ms_s_fwd,stream);
	} // end sim loop
	  
	// clean everything up
      }
      catch (RVLException & e) {
	....
      }
    }
  };

  IWaveAdjSim: public Algorithm {

  private:
    IWaveTreeState iwstate;                       // wraps vector of IWaveStates
    GFD_MODEL gfdm;
    grid2file * rsa;                              // grid output - adjoint of grid input
    file2trace * m_s_adj;                         // trace input - adjoint to trace output
    trace2file * m_s_fwd;                         // trace output
    IWaveCPSim * cps;                             // reference field simulator

  public:
    IWaveAdjSim(PARARRAY & _pars, FILE * _stream, GFD_INIT_FUN minit, 
		unsigned int nsteps, unsigned int nsnaps, unsigned int order)
      : iwstate(_pars,_stream,minit,order),
      pars(_pars), stream(_stream) {
      if (order<1) {
	RVLException e;
	e<<"Error: IWaveAdjSim constructor\n";
	e<<"  order < 1 not permitted\n";
	throw e;
      }
      cps = new IWaveCPSim(_pars,_stream,minit,*(iwstate.ref),nteps,nsnaps);
      minit(&gfdm);
      rsa = new grid2file;
      // need to think up constructor for this thing - combine possible shift
      // to dual grid with rsfread.
      // - use isdyn and specs->fd_set_grid_type to param file2grid
      m_s_adj = new file2trace;
      if (file2trace_construct(m_s_adj,&pars,hdrkey,stream)) {
	// throw exception
      }
      // probably should be a movie option
      // ref simulator
      
    }

    void run() {
      try {

	while (ms_s_adj->xrec <= ms_s_adj->last) {
	  // zero everything
	  iwstate.zero();

	  // initialize sampler - always last component sampled
	  std::vector<IWaveState *> states = iwstate.getStateArray();
	  file2trace_init(ms_s_adj,
			  &(states[states.size()-1]->model),
			  &(iwstate.getPARARRAY()),iwstate.getFILE());
	    
	  // time loop - note that this requres storagee of the sampler
	  // start and stop times in TIMESTEPINDEX form
	  iwstate.setTime(ms_s_adj.stop);
	  TIMESTEPINDEX tsi;
	  // cps gets initial time from source sampler, generated internally
	  cps->setFinalTime(ms_s_adj.stop);
	  while (less_than(ms_s_adj.start,iwstate.getTime())) {
	    // external copy of time stamp
	    tsi_copy(&tsi,&iwstate.getTime());
	    // decrement external copy
	    prev_step(tsi);
	    // set reference time
	    cps->setTargetTime(tsi);
	    // run reference field to reference time
	    cps->run();
	    // trace input
	    file2trace_run(ms_s_adj,&(states[states.size()-1]->model));
	    // backward time step
	    tsa(states);
	    // accumulate adjoint
	    grid2file_run(rsa,states);
	  } //end time loop 
	  
	  // write ouput - note that d, og inputs should be struct members
	  // and initialized by the _init function
	  grid2file_flush(rsa,stream);
	} // end sim loop
	  
	// clean everything up
      }
      catch (RVLException & e) {
	....
      }
    }
  };
    
  /////////////////////// Test Structures ///////////////////////
  
  typedef struct s_TestState0Data {
    int nlvl;
  } TestState0Data;
  
  class TestState0 {
  private:
    int nlvl;
    mutable IWaveTime ts;
    mutable TIMESTEPINDEX lts;
  public:
    TestState0(TestState0Data const & d): nlvl(d.nlvl), ts(lts) {lts.it=0; lts.iv=0; lts.dt=1.0; lts.niv=nlvl;}
    TestState0(TestState0 const & s): nlvl(s.nlvl), ts(s.ts) {}
    void setTime(Time const & t) {
      try {
	ts = t;
      }
      catch (RVLException & e) {
	e<<"\ncalled from TestState0::setTime\n";
	throw e;
      }
    }
    Time & getTime() { return ts; }
    Time const & getTime() const { return ts; }

    void next_step() { 
      TIMESTEPINDEX & dts = ts.getCstruct();
      dts.iv++;
      if (dts.iv >= nlvl) {
	dts.iv=0;
	dts.it++;
      }
    }
    ostream & write(ostream & str) const {
      str<<"TestState0: IWaveTime = \n";
      ts.write(str);
      return str;
    }
  };

  typedef IWaveTreeState<TestState0, TestState0Data> TestTreeState0;

  void TestState0Fun(std::vector<TestState0 *> w);

  
}
