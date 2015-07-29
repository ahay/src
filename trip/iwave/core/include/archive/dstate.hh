#include "except.hh"
#include <deque>
#include "sim.hh"
#include "iwavetime.hh"

namespace TSOpt {

  // array of IWaveStates with ref object of same
  // type pointing to first half of array
  class IWaveTreeState {
  private:
    bool own;
    std::vector<IWAVE *> sa;
    std::vector<RDOM *> rd;
    IWaveTreeState * ref;
    IWaveTreeState();
    IWaveTreeState(IWaveTreeState const &);
    IWaveTreeState(std::vector<IWAVE *> sv)
      : own(false), sa(iwave_max(sv.size()/2,0)),
	rd(iwave_max(sv.size()/2,0)) {
      if (sv.size()>1) 
	for (int i=0;i<sv.size()/2;i++) {
	  sa[i]=sv[i];
	  rd[i]=(sv[i]->model).ld_c;
	}
      ref=NULL;
    }
  public:
    IWaveTreeState(PARARRAY & _pars, FILE * _stream,
		   GFDM_INIT_FUN _minit, int n=0): own(true), sa(pow2(n)) {
      try {
	for (size_t i=0; i<sa.size(); i++) {
	  sa[i]=new IWAVE(_pars,_stream,_minit);
	  rd[i]=(sa[i]->model).ld_c;
	}
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
    
    virtual std::vector<IMODEL *> & getStateArray() { return sa; }
    virtual std::vector<IMODEL *> const & getStateArray() const { return sa; }

    ostream & write(ostream & str) const {
      str<<"IWaveTreeState, length "<<sa.size()<<"\n";
      return str;
    }
  };

    
  /* assumption: each sampler (explicit for traces, implicit for grids) 
     samples to/from a unique RARRAY in the IMODEL RDOM, and is tied to 
     a unique file built on a unique prototype 
  */
  class FileToGrid {

  private:

    std::vector<TASK_RELN *> tr;
    
    int nrec;  // total number of records in simulation
    int first; // index of first record for this process
    int last;  // index of last record for this process

    PARARRAY const & par;
    FILE * stream;

  public:

    FileToGrid(std::vector<TASK_RELN *> _tr,
	       PARARRAY const & _par,
	       FILE * _stream) 
      : tr(_tr.size()),
	par(_par), nrec(0), irec(0), xrec(0), samplers(_indices.size(), NULL),
	stream(_stream) {
      try {

	/* deep copy of TASK_RELN array */
	for (int i=0;i<tr.size();tr++) {
	  tr[i]=new TASK_RELN;
	  deep_copy(tr[i],_tr[i]);
	}

	/* sanity */

	int init = 0;
	for (int n = 0; n< tr.size(); n++) {
	  
	  /* postcondition of call to sampler_construct: tracegeom object
	     has all memory allocated, and basic data set - esp. 
	     nrec - number of records
	     irec - current record number
	     xrec - next record nubmer
	     ntr[it] - number of traces in record it
	  */

	  SAMPLER * trace = new SAMPLER;
	  // hack to deal with present needs of sampler
	  // should hack sampler to eliminate these
	  IPNT sindex;
	  RPNT mindex;
	  RPNT scoord;
	  for (int i=0; i< RARR_MAX_NDIM; i++) {
	    sindex[i]=0;
	    mindex[i]=1.0f;
	    scoord[i]=0.0f;
	  }
	  sindex[0]=indices[n];
	  // if sampler construction fails, don't increment count
	  // and don't add pointer to list
	  if (sampler_construct(&trace,
				&pars,
				sindex,
				mindex,
				scoord,
				input[n],
				hdrkey[n],
				datakey[n],
				stream)) {
	    fprintf(stream,"NOTE: FileToGrid constructor\n");
	    fprintf(stream,"      no trace sampler construction on array %d\n",indices[n]);
	    delete trace;
	  }
	  else {
	    sampler[n]=trace;
	    init++;
	  
	    /* more sanity check - besure that the number of records 
	       is the same for all trace samplers */
	    if (nrec) {
	      if ((nrec  == trace.t.tg.nrec) &&
		  (first == trace.t.tg.first) &&
		  (last  == trace.t.tg.last)) nrec++;
	      else {
		RVLException e;
		e<<"Error: trace sampler on input "<<n<<" has different number\n";
		e<<"  of records from first trace sampler encountered,\n";
		e<<"  or first or last record indices for this process \n";
		e<<"  are inconsistent with those extraced from first \n";
		e<<"  trace sampler\n";
		throw e;
	      }
	    }
	    else {
	      // first trace sampler
	      nrec  = trace.t.tg.nrec;
	      first = trace.t.tg.first;
	      last  = trace.t.tg.last;
	    }
	  }
	}
	
	/* more sanity: need at least one trace sampler if there is to
	 be more than one record. Reason: the roles of trace and
	 grid sampling are slightly asymmetric, because the
	 simulation may have many output records (shots) but only a
	 single set of input coefficient fields. In that case the
	 number of records to simulate MUST come from the trace
	 files, as there is no info about it anywhere else in the
	 simulation. If we come to output non-trace data exclusively
	 in some problem, then it will be necessary to figure it
	 out. Notice that even if there is only a source file, it
	 takes the form of trace data and is sufficient to specify
	 the number of records.*/ 
	if (init < 1) {
	  RVLException e;
	  e<<"Error: FileToGrid\n";
	  e<<"  require at least one trace sampler to define record structure\n";
	  e<<"  example: an SU file (SEGY traces) suffices\n";
	  throw e;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from FileToGrid constructor\n";
	throw e;
      }
    }

    void load(IMODEL & m, int panelindex) {
      try {
	for (int n=0; n<indices.size(); n++) {
	  if (input[n]) {

	    // first test rsf case

	    /*==================================================
	      ================= DECLARATIONS ===================
	      ==================================================*/
	    
	    /* tmp storage for read */
	    RARR Btmp;
	    
	    /* other workspace */
	    int i,j,m; /* counters */
	    int dim;         /* problem dimension */
	    IPNT ran;        /* axis lens */
	    IPNT rags;       /* axis starts for read */
	    IPNT rage;       /* axis ends for read */
	    IPNT raVgs;      /* axis starts for output */
	    IPNT raVge;      /* axis ends for output */
	    unsigned long ntot;        /* total num of words in phys dom */
	      
	    /* error flag */
	    int err=0;
	    IPNT ip;
	    ireal q;
	      
	    /*==================================================
	      =============== END DECLARATIONS =================
	      ==================================================*/  
	      
	    /* axis lens, starts - allocated, not virtual */
	    rd_ndim(&dom,indices[n],&dim);
	    rd_size(&dom,indices[n],ran);
	    rd_gse(&dom, indieces[n],raVgs,raVge);
	    ntot=1;
	    for (int i=0;i<dim;i++) ntot*=ran[i];
	    // zero non-physical dims
	    for (int i=dim; i<RARR_MAX_NDIM; i++) {
	      raVgs[i]=0;
	      raVge[i]=0;
	      rags[i]=0;
	      rage[i]=0;
	      ran[i]=1;
	    }
	    // zero target rarr
	    ra_a_zero(&(m.ld_a._s[indices[n]]));

	    /* extract grid type */
	    IPNT gtype[RDOM_MAX_NARR];
	    if (((FD_MODEL *)(m.specs))->set_grid_type(stream,dim,gtype)) {
	      RVLException e;
	      e<<"Error: FileToGrid::load from FD_MODEL::set_grid_type\n";
	      throw e;
	    }

	    /* bump up array limits on shifted dimensions */
	    for (int i=0; i<dim; i++) {
	      rags[i]=raVgs[i];
	      if (gtype[indices[n]][i]) {
		rage[i]=raVge[i]+1;
		ran[i]++;
	      }
	      else rage[i]=raVge[i];
	    }
	      
	    /* create read buffer */
	    if (ra_create(&Btmp,dim,rags,rage)) {
	      RVLException e;
	      e<<"Error: FileToGrid::load from ra_create\n";
	      e<<"  failed to create temp buffer\n";
	      throw e;
	    }

	    /* suck the data in */
	    if (rsfread(Btmp._s0,rags,ran,dfile,1,stream,panelindex)) {
	      ra_destroy(&Btmp);
	      fprintf(stream,"NOTE: FileToGrid::load\n");
	      fprintf(stream,"      no rsf load on index %d\n",indices[n]);
	    }
	    else {

	      /* NOTE: THIS SEGMENT ASSUMES THAT 3 <= RARR_MAX_NDIM <= 5 */

	      /* extract data pointer (allocated) to correct index */
	      float * p = m.ld_a._s[indices[n]]._s0;

	      /* write data onto target rarray, left-justified */
#if RARR_MAX_NDIM > 4
	      for (int i4=raVgs[4]; i4<=raVge[4]; i4++) {
#endif
#if RARR_MAX_NDIM > 3
		for (int i3=raVgs[3]; i3<=raVge[3]; i3++) {
#endif
		  for (int i2=raVgs[2]; i2<=raVge[2]; i2++) {
		    for (int i1=raVgs[1]; i1<=raVge[1]; i1++) {
		      for (int i0=raVgs[0]; i0<=raVge[0]; i0++) {
			// sum over vertices of squashed cube
#if RARR_MAX_NDIM > 4			
			for (ishift[4]=0;ishift[4]<=gtype[n][4];ishift[4]++) {
#endif
#if RARR_MAX_NDIM > 3			
			  for (ishift[3]=0;ishift[3]<=gtype[n][3];ishift[3]++) {
#endif
			    for (ishift[2]=0;ishift[2]<=gtype[n][2]; ishift[2]++) {
			      for (ishift[2]=0;ishift[2]<=gtype[n][2]; ishift[1]++) {
				for (ishift[2]=0;ishift[2]<=gtype[n][2]; ishift[0]++) {
				  p
#if RARR_MAX_NDIM > 4
				    [i4]
#endif
#if RARR_MAX_NDIM > 3
				    [i3]
#endif
				    [i2][i1][i0]
				    
				    += pow2(-RARR_MAX_NDIM)*Btmp._s0
#if RARR_MAX_NDIM > 4
				    [i4+ishift[4]]
#endif
#if RARR_MAX_NDIM > 3
				    [i3+ishift[3]]
#endif
				    [i2+ishift[2]][i1+ishift[1]][i0+ishift[0]]
				    
				    }
			      }
			    }
#if RARR_MAX_NDIM > 3
			  }
#endif
#if RARR_MAX_NDIM > 4
			}
#endif
		      }
		    }
		  }
#if RARR_MAX_NDIM > 3
		}
#endif
#if RARR_MAX_NDIM > 4
	      }
#endif
	    }

	  
	    // su case
	    if (samplers[n]) {
	      SAMPLER * s = dynamic_cast<SAMPLER *>(samplers[n]);
	      if (s && sampler_init(s,&m,&par,stream)) {
		if (s->t.tg.irec != panelindex) {
		  RVLException e;
		  e<<"Error: FileToGrid::load\n";
		e<<"  sampler for SU datakey = "<<datakey[n]<<" not synched with sim\n";
		e<<"  sampler record number = "<< s->t.tg.irec<<"\n";
		e<<"  sim record number     = "<< panelindex <<"\n";
		throw e;
	      }
	    }
	    else {
	      RVLException e;
	      e<<"Error: FileToGrid::load\n";
	      e<<"  either sampler for SU datakey = "<<datakey[n]<<" is either\n";
	      e<<"  not constructed or not initialized\n";
	      throw e;
	    }
	  }
	    
	  else {
	    RVLException e;
	    e<<"Error: FileToGrid constructor\n";
	    e<<"proto file "<<hdrkey[n]<<" of unidentified type\n";
	    throw e;
	  }
	}	      
      }
      catch (RVLException & e) {
	e<<"\ncalled from FileToGrid::load\n";
	throw e;
      }
    }

    void sample(IMODEL & m) {
      // no-op for non-trace. This could be changed: grid data could be sampled
      // at initial time. However for now leave in load function.
      
      try {

	for (int n=0; n<indices.size(); n++) {
	  if (samplers[n])) {
	    SAMPLER * s = NULL;
	    s = dynamic_cast<SAMPLER *>(samplers[n]);
	    if (!(s && sampler_run(s,&m))) {
	      RVLException e;
	      e<<"Error: FileToGrid::sample\n";
	      e<<"  either sampler not constructed, or sampler_run failed\n";
	      e<<"  index = "<<n<<" data key = "<<datakey[n]<<"\n";
	      throw e;
	    }
	  }
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from FileToGrid::sample\n";
	throw e;
      }
    }
    
    void save(IMODEL & m) {
    }
    
    ostream & write(ostream & str) const {
      
    }

  };

  /* input IMODELS for fwd - output for adj: 0, 2^0, 2^1,..2^{n-1}, output for 2^{n-1} 
     for basic modeling (0th deriv): create FileToGrid for sa[0] - input and output
     for 1st deriv: for sa[0] (input), sa[1] (input and output)
     for 2nd deriv: for sa[0] (input), sa[1] (input), sa[2] (input), and sa[3] (output)
     for 3rd deriv: for 

   */

  class IWaveSim: public Algorithm {

  private:
    IWaveTreeState iwstate;                       // wraps vector of IMODELs
    GFD_MODEL gfdm;                               // forward time step function
    FileToGrid * fg;                              // i/o
    PARARRAY & pars;
    FILE * stream;

  public:
    IWaveSim(int order, PARARRAY & _pars, FILE * _stream, GFD_INIT_FUN minit)
      : iwstate(_pars,_stream,minit,order), pars(_pars), stream(_stream) {
      minit(&gfdm);
      fg = new FileToGrid(gfdm.getIndices(),
			  gfdm.getHdrkey(),
			  gfdm.getDatakey(),
			  PARARRAY const & _par,
			FILE * _stream) 
      // probably should be a movie option
      }

  void run() {
    try {
      // simulation loop
      while (ms_s_fwd->xrec <= ms_s_fwd->last) {
	// zero everything
	iwstate.zero();
	// load sim data - used to be iwave_static_init
	fg.load();
	// initialize sampler - always last component sampled
	std::vector<IMODEL *> states = iwstate.getStateArray();
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
  IWaveTreeState iwstate;                       // wraps vector of IMODELs
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
	std::vector<IMODEL *> states = iwstate.getStateArray();
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
	} //end time loop 
	// write ouput - note that d, og inputs should be struct members
	// and initialized by the _init function
	grid2file_run(rsa,states,stream);
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
