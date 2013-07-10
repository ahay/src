#ifndef __IWAVE_OP
#define __IWAVE_OP

#include "alg.hh"
#include "op.hh"
#include "ocdc.hh"
#include "CPsim.hh"
#include "iwstack.hh"
#include "parserdecl.hh"

namespace TSOpt {

  //using namespace RVLAlg;
  using RVLAlg::ListAlg;
  using RVL::Space;
  using RVL::ProductSpace;
  using RVL::ConstContainer;
  using RVL::STRING_PAIR;
  using RVL::Vector;
  using RVL::FunctionObject;
  using RVL::Operator;
  using RVL::Writeable;
  using RVL::AssignParams;
  
  /** SamplePolicy creates Sampler, incorporating all input and output
      grid sampling; SimPolicy creates Sim type

      One each for fwd, lin, and adj modes - each samples, steps the
      corresponding field. For Lin and Adj modes, also need separate
      choices for stepping reference field - these are LinFwdSimPolicy
      and AdjFwdSimPolicy
  */
  template
  <
    class FwdSamplerPolicy,
    class LinSamplerPolicy,
    class AdjSamplerPolicy,
    class FwdSimPolicy,
    class LinFwdSimPolicy,
    class LinSimPolicy,
    class AdjFwdSimPolicy,
    class AdjSimPolicy
  >
  class IWaveOp: public Operator<ireal>, 
		 public FwdSamplerPolicy,
		 public LinSamplerPolicy,
		 public AdjSamplerPolicy,
		 public FwdSimPolicy,
		 public LinFwdSimPolicy,
		 public LinSimPolicy,
		 public AdjFwdSimPolicy,
		 public AdjSimPolicy
  {
      
  private:
    
    mutable FILE * stream;              /* output stream            */
    // Change 11.05.10 - copy mutable, const ref copy
    mutable PARARRAY * lpars;           /* parameter array workspace */
    // change 23.06.10 - ref mutable too, to facilitate 
    // under-the-table definiton of operator on source space
    mutable PARARRAY * pars;            /* parameter array ref copy */

    mutable GFDM_INIT_FUN minit;

    Space<ireal> const & dom;
    Space<ireal> const & rng;

    // verbosity control
    int dump_steps;
    int dump_pars;
    int dump_term;

    // other verbosity control handled within iwave code
    
    IWaveOp();
      
  protected:
      
    void apply(const Vector<ireal> & x, 
	       Vector<ireal> & y) const {

      try {

	SpaceTest(this->getDomain(),x,"TSOpt::IWaveOp::apply (dom)");
	SpaceTest(this->getRange(),y,"TSOpt::IWaveOp::apply (rng)");

	// Change 11.05.10 - lpars, not pars, referenced in IWaveLinState
	// copy over reference (pars), add parameters from input, output
	// data 
	ps_copy(&lpars,*pars);
	AssignParams ap(*lpars,stream);
	y.eval(ap);
	x.eval(ap);

	// Change 11.05.10 - print par file used in sim
	if (dump_pars) {
	  fprintf(stream,"PARAM ARRAY CREATED IN IWOP::APPLY\n");
	  ps_printall(*lpars,stream);
	  fflush(stream);
	}

	/* zero output */
	y.zero();

	if (dump_steps) {
	  fprintf(stream,"IWOP::APPLY -> STEP\n");
	  fflush(stream);
	}

	IWaveState iw(*lpars,stream,minit);

	IWaveStep iwstep(iw);

	if (dump_steps) {
	  fprintf(stream,"IWOP::APPLY -> SAMPLER\n");
	  fflush(stream);
	}

	/* create sampler */
	/* sampler needs to own
	   - method init() returning bool - true if another sim is ready
	   - algs pre(), post(), and term(), the latter being a time term
	   - start() should return the TSIndex rep of start time
	   - term() should include any final output with true return
	   - getCoords(RPNT & src) method to update source coords (implement
	     Coords interface as well as Sampler interface)
	*/
	Sampler * s = NULL;
	if (!(s=FwdSamplerPolicy::create(iw))) {
	  RVLException e;
	  e<<"Error: IWaveOp::apply\n";
	  e<<"failed to create sampler from policy\n";
	  throw e;
	}
	if (dump_steps) {
	  fprintf(stream,"IWOP::APPLY -> INIT\n");
	  fflush(stream);
	}

	/* Note: IWaveStaticInit object stores reference to sampler,
	   for cases in which static fields depend on output record
	   number (surface-oriented extension)
	*/
	IWaveStaticInit iwstatinit(iw,*s);
	/* initializes dynamic fields */
	IWaveDynamicInit iwdyninit(iw);
	/* advances time step data structure */
	IWaveNextTime<IWaveState> iwnexttime(iw);

	if (dump_steps) {
	  fprintf(stream,"IWOP::APPLY -> STATINIT\n");
	  fflush(stream);
	}

	iwstatinit.run();
	
	if (dump_steps) {
	  fprintf(stream,"IWOP::APPLY -> LOOP\n");
	  fflush(stream);
	}

	/* sim loop */
	while (s->init()) {	 

	  // this call may be a fake - must use the Sampler 
	  // interface, but apps may write to the FILE * stream

	  if (dump_term) {
	    s->write(cout);
	  }

	  // list includes 
	  // s->pre: source injection
	  // s->post: trace/movie projection
	  // mod of 10.04.12: first inject source, then update 
	  // time step, then sample
	  ListAlg iwpost1(s->pre(),iwnexttime);
	  // append time advance to list
	  ListAlg iwpost(iwpost1,s->post());
	  // prepend time step to create final list alg for loop
	  TimeStepList<IWaveState> iwlist(iwstep,iwpost,false);
	 
	  // set start and end times for simulation - transfers information
	  // from acquisition geometry to state
	  iwdyninit.setTimerange(s->start(),s->end());	
	  // initialize dynamic fields  
	  iwdyninit.run(); /* must run it before build iwstack */
	  
	  // for convenience, encapsulate sim data in a standard object
	  StdSimData<IWaveState> sdata(iwlist,s->term(),iwdyninit);

	  // use policy to create simulator; run; delete simulator 
	  // object; flush results to archival storage
      	  Sim<IWaveState> * ts = NULL;
	  if (!(ts=FwdSimPolicy::create(sdata))) {
	    RVLException e;
	    e<<"Error: IWaveOp::apply\n";
	    e<<"failed to create sim from policy\n";
	    throw e;
	  }	    
	  ts->run();
	  delete ts;
	  s->flush();

	  // for surface-oriented extension, reinitialize static
	  // fields for (possible) next simulation - no-op if 
	  // model is not extended, or if record number is out of range
	  iwstatinit.run();
	}
	delete s;

	if (dump_steps) {
	  fprintf(stream,"IWOP::APPLY -> EXIT\n");
	  fflush(stream);
	}

	// placed to force all output to be written to disk 
	// before any is copied
#ifdef IWAVE_USE_MPI
	MPI_Barrier(retrieveGlobalComm());
#endif
      }
      catch (RVLException & e) {
	e<<"\ncalled from TSOpt::IWaveOp::apply\n";
	throw e;
      }
    }
      
    void applyDeriv(const Vector<ireal> & x, 
		    const Vector<ireal> & dx,
		    Vector<ireal> & dy) const {

      try {

	// sanity test arguments
	SpaceTest(this->getDomain(),x,"TSOpt::IWaveOp::applyDeriv (dom, ref)");
	SpaceTest(this->getDomain(),dx,"TSOpt::IWaveOp::applyDeriv (dom, pert)");
	SpaceTest(this->getRange(),dy,"TSOpt::IWaveOp::applyDeriv (rng)");

	// Change 11.05.10 - lpars, not pars, referenced in IWaveLinState
	// copy over reference (pars), add parameters from input, output
	// data 
	ps_copy(&lpars,*pars);

	/* write filenames to internal table: reference fields */
	AssignParams ap(*lpars,stream);
	x.eval(ap);

	/* write filenames to internal table: pert fields */
	string dstr="d"; 
	string nstr="";
	AssignParams dap(*lpars,dstr,nstr,stream); 
	dx.eval(dap);

	/* write filenames to internal table: output fields */
	string lstr="l"; 
	AssignParams lap(*lpars,lstr,nstr,stream); 
	dy.eval(lap);

	if (dump_pars) {
	  fprintf(stream,"PARAM ARRAY CREATED IN IWOP::APPLYDERIV\n");
	  ps_printall(*lpars,stream);
	  fflush(stream);
	}

	/* zero output */
	dy.zero();

	if (dump_steps) {
	  fprintf(stream,"IWOP::APPLYDER -> STEP\n");
	  fflush(stream);
	}

	/* tangent state */
	/* steps are actually the same, but applied to different
	   dynamic fields */
	IWaveLinState iw(*lpars,stream,minit);
	IWaveStep iwfwdstep(iw);
	IWaveLinStep iwlinstep(iw);

	if (dump_steps) {
	  fprintf(stream,"IWOP::APPLYDER -> SAMPLER\n");
	  fflush(stream);
	}

	/* sampler needs to own
	   - method init() returning bool - true if another sim is ready
	   - algs pre(), post(), and term(), the latter being a time term
	   - algs linpre(), linpost()
	   - start() should return the TSIndex rep of start time
	   - term() should include any final output with true return
	*/
	LinSampler * s = NULL;
	if (!(s=LinSamplerPolicy::create(iw))) {
	  RVLException e;
	  e<<"Error: IWaveOp::applyDer\n";
	  e<<"failed to create Sampler for linearized state from policy\n";
	  throw e;
	}

	if (dump_steps) {
	  fprintf(stream,"IWOP::APPLYDER -> INIT\n");
	  fflush(stream);
	}

	/* similar for inits - as for steps, mediated by getIWAVE */
	IWaveStaticInit iwstatinit(iw,*s);
	IWaveDynamicInit iwdyninit(iw);
	//	IWaveStaticInit<IWaveLinState> iwlinstatinit(iw,*c);
	IWaveLinDynamicInit iwlindyninit(iw);
	IWaveNextTime<IWaveState> iwnexttime(iw);
	IWaveNextTime<IWaveLinState> iwlinnexttime(iw);

	if (dump_steps) {
	  fprintf(stream,"IWOP::APPLYDER -> REF-STATE STATINIT\n");
	  fflush(stream);
	}

	iwstatinit.run();

	if (dump_steps) {
	  fprintf(stream,"IWOP::APPLYDER -> LOOP\n");
	  fflush(stream);
	}

	/* sim loop */
	while (s->init()) {

	  if (dump_term) {
	    s->write(cout);
	  }

	  /* action lists for reference, linearized timesteps */
	  // mod of 10.04.12: add source first, then update time step
	  ListAlg iwpost(s->pre(),iwnexttime);
	  // big change 17.01.12: replace s->linpre (call to acoustic born source code) 
	  // with generic IWaveDModStep 
	  // mod of 10.04.12: update time step first, then sample
	  ListAlg iwlinpost(iwlinnexttime,s->linpost());

	  TimeStepList<IWaveState> iwfwdlist(iwfwdstep,iwpost,false);

	  TimeStepList<IWaveLinState> iwlinlist(iwlinstep,iwlinpost,false);

	  /* set start times */
	  iwdyninit.setstartTime(s->start());
	  iwlindyninit.setTime(s->start());

	  /* run iwdyninit to assure init of ref field */
	  iwdyninit.run(); /* must run it before build iwstack */
	 
	  /* build fwd simulator - NOTE init not used SHOULD FIXXXXXX */
	  StdSimData<IWaveState> sdata(iwfwdlist,s->term(),iwdyninit);	 
	  Sim<IWaveState> * tsfwd = NULL;

	  if (!(tsfwd=LinFwdSimPolicy::create(sdata))) {
	    RVLException e;
	    e<<"Error: IWaveOp::applyDer\n";
	    e<<"failed to create sim from policy\n";
	    throw e;
	  }	    

	  /* couple fwd sim, lin step */
	  //	  TargetCurrTimeStep<IWaveLinState,IWaveState> synch(*tsfwd,iwlinstep);
	  IWaveTargetRefTime linsynch(*tsfwd,iwlinstep,true);

	  /* build action list incorporating fwd update */
	  TimeStepList<IWaveLinState> linlist(iwlinlist,linsynch,true);

	  /* build lin sim */
	  StdSimData<IWaveLinState> slindata(linlist,s->linterm(),iwlindyninit);
	  Sim<IWaveLinState> * tslin = NULL;
	  if (!(tslin=LinSimPolicy::create(slindata))) {
	    RVLException e;
	    e<<"Error: IWaveOp::applyDer\n";
	    e<<"failed to create lin sim from policy\n";
	    throw e;
	  }	    
	  /* run */
	  tslin->run();
	  /* clean up */
	  delete tslin;
	  delete tsfwd;
	  /* write data */
	  s->flush();

	  iwstatinit.run();
	}
	delete s;

	if (dump_steps) {
	  fprintf(stream,"IWOP::APPLYDER -> EXIT\n");
	  fflush(stream);
	}

	// placed to force all output to be written to disk 
	// before any is copied
#ifdef IWAVE_USE_MPI
	MPI_Barrier(retrieveGlobalComm());
#endif
      }

      catch (RVLException & e) {
	e<<"\ncalled from TSOpt::IWaveOp::applyDeriv\n";
	throw e;
      }

    }
      
    void applyAdjDeriv(const Vector<ireal> & x, 
		       const Vector<ireal> & dy,
		       Vector<ireal> & dx) const {

      try{

	if (dump_steps) {
	  fprintf(stream,"IWOP::APPLYADJ -> START\n");
	  fflush(stream);
	}	

	SpaceTest(this->getDomain(),x,"TSOpt::IWaveOp::applyAdjDeriv (dom, ref)");
	SpaceTest(this->getRange(),dy,"TSOpt::IWaveOp::applyAdjDeriv (rng)");
	SpaceTest(this->getDomain(),dx,"TSOpt::IWaveOp::applyAdjDeriv (dom, pert)");
	
	// Change 11.05.10 - lpars, not pars, referenced in IWaveLinState
	// copy over reference (pars), add parameters from input, output
	// data 
	ps_copy(&lpars,*pars);

	/* write filenames to internal table*/
	AssignParams ap(*lpars,stream);

	/* reference (forward) fields */
	x.eval(ap);
	
	/* input data residual*/
	string rstr="r";
	string nstr="";
	AssignParams rap(*lpars,rstr,nstr,stream);
	dy.eval(rap);

        /* output adj pert model*/
	string mstr="m";
	AssignParams map(*lpars,mstr,nstr,stream);
	dx.eval(map);

	if (dump_pars) {
	  fprintf(stream,"PARAM ARRAY CREATED IN IWOP::APPLYADJDERIV\n");
	  ps_printall(*lpars,stream);

	  fprintf(stream,"DEBUG OUTPUT: FILE MGR STATUS\n");
	  iwave_fprintall(stream);
	  fflush(stream);
	}

	/* zero output */
	dx.zero();

	if (dump_steps) {
	  fprintf(stream,"IWOP::APPLYADJ -> STEP\n");
	  fflush(stream);
	}

       	/* tangent state */
	
 	/*change time-step info for backward field (i.e.,
	  iw.IWaveLinState::getIWAVE()), such that dt_bwd = - dt_fwd,
	  lam_fwd = -lam_bwd.  

          ------ NOTE: Currently the above is done in
          ASGAdjSampler::init() ------
	*/
	IWaveLinState iw(*lpars,stream,minit);
	IWaveStep iwfwdstep(iw);
 	IWaveLinStep iwbwdstep(iw,false);

	/* false indicates bwd-step */

	if (dump_steps) {
	  fprintf(stream,"IWOP::APPLYADJ -> SAMPLER\n");
	  fflush(stream);	
	}
 	/* create sampler for adjoint state computation */
 	/* adj sampler needs to own 

	   - method init() returning bool - true if another sim is
	   ready

	   - algs pre(), post(), and term(), components of building
	   fwd sim object

	   - algs linpre(), linpost(), and linterm(), components for
	   building adj computation linpre(): imaging condition,
	   linpost(): inserting src for bwd sim, linterm(): terminate
	   bwd loop (i.e., adj computation)

	   [BIG CHANGE 17.01.12: use iwdmodstep instead of s->linpre -
	   decouple Born source / imaging condition from sampler]

	   - start() should return the TSIndex rep of start time

	   - end() should return the TSIndex rep of final time (final
	   recording time)
 	*/
 	LinSampler * s = NULL;

 	if (!(s=AdjSamplerPolicy::create(iw))) {
 	  RVLException e;
 	  e<<"Error: IWaveOp::applyAdjDeriv\n";
 	  e<<"failed to create Sampler for IWaveLinstate from policy\n";
 	  throw e;
 	}

	/* extract number of snaps for checkpointing - default = 1 */
	int nsnaps = 3;
	ps_flint(map.get(),"nsnaps",&nsnaps);
        /* extract maximum number of stack buffers - default = 0 */
        int max_num_buffer = 0;
        ps_flint(*pars,"maxbuffernum",&max_num_buffer);

	if (dump_steps) {
	  fprintf(stream,"IWOP::APPLYADJ -> INIT\n");
	  fflush(stream);	
	}

 	/* static & dynamic initializing algorithms */
 	IWaveStaticInit iwfwdstatinit(iw,*s);
 	IWaveDynamicInit iwfwddyninit(iw);
 	IWaveLinDynamicInit iwbwddyninit(iw);

	IWaveNextTime<IWaveState> iwfwdnexttime(iw);
	IWaveNextTime<IWaveLinState> iwbwdnexttime(iw,false);
	
	if (dump_steps) {
	  fprintf(stream,"IWOP::APPLYADJ -> FWD-STATINIT\n");
	  fflush(stream);
	}	
	iwfwdstatinit.run();

	if (dump_steps) {
	  fprintf(stream,"IWOP::APPLYADJ -> LOOP\n");
	  fflush(stream);	
	}
 	/* sim loop: (outer loop through acquisition parameter, such
	   as source coords) */
 	while (s->init()) {

	  if (dump_term) {
	    s->write(cout);
	  }

 	  /* action lists for forward, backward timesteps */
	  ListAlg iwpost(s->pre(),iwfwdnexttime);
      	  TimeStepList<IWaveState> iwfwdlist(iwfwdstep,iwpost,false);

	  //////////////////////////////////////////////////////////////
	  // notes from ASG_AdjSampler:
	  // linpre = adjtrs = adjoint trace sampler
	  // linpost = abs = no-alg (now!)
	  // 17.01.12: replace s->linpost with iwdmodstep
	  //////////////////////////////////////////////////////////////
	  //
	  // derivation: lin fwd loop is schematically
	  //
	  // time = start;
	  // while time<stop
	  //   synch with reference state
	  //   update state:
	  //     call update function
	  //     time++
          //   sample state
	  // end
	  // 
	  // so lin adj loop is 
	  // 
	  // time = stop;
	  // while time > start  
	  //   synch with reference state
	  //   adjoint sample state
	  //   update state:
          //     time--
          //     call update function 
	  // end
	  // 
	  // The components are
	  // synch with reference state: adjsynch [couples forward "restore" 
          //    simulator to adjoint step]
	  // state update: iwbwdstep
	  // adjoint sampler: s->linpre();
          // time decrement operator: iwbwdnexttime
	  //
	  //////////////////////////////////////////////////////////////////

	  // adjoint sample, then decrement time
	  ListAlg iwbwdpre(s->linpre(),iwbwdnexttime);
	  // insert previous before adj step
	  TimeStepList<IWaveLinState> iwbwdlist(iwbwdstep,iwbwdpre,true);
	 
 	  /* set starting time for fwd step*/
 	  iwfwddyninit.setstartTime(s->start());
	  /* set starting time for bwd step */
	  iwbwddyninit.setTime(s->end());

	 
	  iwfwddyninit.run(); /* must run it before build iwstack */
  
	  int tmp_its, tmp_ite;
	  tmp_its = tmp_ite = 0;
	  try {
	    TSIndex const & tmp_tstart = dynamic_cast<TSIndex const &>(s->start());
	    TSIndex const & tmp_tend = dynamic_cast<TSIndex const &>(s->end());

	    tmp_its = tmp_tstart.getCstruct().it;
	    tmp_ite = tmp_tend.getCstruct().it;
    	  }
	  catch (const std::bad_cast& e) {
	    std::cerr << e.what() << std::endl;
	    std::cerr << "IWaveOp::applyAdjDeriv, the object returned by "
		      <<"s->start() or s->end() is not of type TSIndex" 
		      << std::endl;
	  }

	  IWaveTimeStack * iwstack = new IWaveTimeStack(iw,tmp_its,max_num_buffer);
	  /* build fwd simulator */
	  StdSimData<IWaveState> 
	    sdata(iwfwdlist,s->term(),iwfwddyninit,iwstack,
		  tmp_ite-tmp_its + 1,nsnaps);	 
       	  Sim<IWaveState> * tsfwd = NULL;
 	  if (!(tsfwd=AdjFwdSimPolicy::create(sdata))) {
 	    RVLException e;
 	    e<<"Error: IWaveOp::applyAdjDeriv\n";
 	    e<<"failed to create sim from policy\n";
 	    throw e;
 	  }	   

	  // synchronizing object: ref sim -> adj step
	  IWaveTargetRefTime adjsynch(*tsfwd,iwbwdstep,false);

	  // synchronization precedes other backward steps
 	  TimeStepList<IWaveLinState> adjlist(iwbwdlist,adjsynch,true);
	  
 	  // put together action list, terminator, and initializer to 
	  // make simulator
 	  StdSimData<IWaveLinState> sadjdata(adjlist,s->linterm(),iwbwddyninit);
 	  Sim<IWaveLinState> * tsadj = NULL;
	  if (!(tsadj=LinSimPolicy::create(sadjdata))) { 	 
	    RVLException e;
 	    e<<"Error: IWaveOp::applyAdjDer\n";
 	    e<<"failed to create adj sim from policy\n";
 	    throw e;
 	  }	  
  
	  //	  tsadj->write(cerr);
 	  /* run */
 	  tsadj->run();

 	  /* clean up */
 	  delete tsadj;
	  delete tsfwd;
	  if (iwstack) delete iwstack;

	  /* for extended models: write adjoint output for this panel,
	     reinitialize static arrays for next panel. if model is not extended, 
	     both of these are no-ops. Barrier required for synthesis */
	  s->flush();
	  iwfwdstatinit.run();

 	}

	/* non-extended modeling: flush at the end of the simulation
	   use barrier to ensure that all groups have finished
	   accumulating partial adjoint output before it is reduced to
	   global adjoint output
	*/
	//	s->flush();
 	delete s;

	// placed to force all output to be written to disk 
	// before any is copied
#ifdef IWAVE_USE_MPI
	MPI_Barrier(retrieveGlobalComm());
#endif
	if (dump_steps) {
	  fprintf(stream,"IWOP::APPLYADJ -> EXIT\n");
	  fflush(stream);	
	}
      }
      catch (RVLException & e) {
	e<<"\n called from TSOpt::IWaveOp::applyAdjDeriv\n";
	throw e;
      }

    }
      
    Operator<ireal> * clone() const { 
      return new IWaveOp<
	FwdSamplerPolicy,
	LinSamplerPolicy,
	AdjSamplerPolicy,
	FwdSimPolicy,
	LinFwdSimPolicy,
	LinSimPolicy,
	AdjFwdSimPolicy,
	AdjSimPolicy
	>(*this);
    }
      
  public:

    // Change 11.05.10: IWaveLinState stores reference to mutable
    // PARARRAY lpars, not to reference pars. lpars updated before use
    // in every apply method, by copying pars.
    IWaveOp(Space<ireal> const & _dom,
	    Space<ireal> const & _rng,
	    PARARRAY _pars, FILE * _stream,
	    GFDM_INIT_FUN _minit)
      : dom(_dom), rng(_rng), 
	stream(_stream),
	pars(NULL),
	lpars(NULL),
	minit(_minit),
        dump_steps(0),
	dump_pars(0),
	dump_term(0)
    {
      /* copy parray, extract filenames from input (which may
	 have multiple components) and output (ditto), presuming 
	 that both are space types which can also be cast to ConstContainers.
      */

      int err=0;

      /* copy input par array to data member */
      if (err=ps_copy(&pars,_pars)) {
	RVLException e;
	e<<"Error: IWOP constructor from ps_copy, err="<<err<<"\n";
	throw e;
      }

      // set dump controls
      ps_flint(*pars,"dump_steps",&dump_steps);
      ps_flint(*pars,"dump_pars",&dump_pars);
      ps_flint(*pars,"dump_term",&dump_term);

      /* see what we've got */
      if (dump_pars) {
	fprintf(stream,"PARS IN IWOP CONSTRUCTOR\n");
	ps_printall(*pars,stream);
      }
    }

    IWaveOp(IWaveOp
	    <
	    FwdSamplerPolicy,
	    LinSamplerPolicy,
	    AdjSamplerPolicy,
	    FwdSimPolicy,
	    LinFwdSimPolicy,
	    LinSimPolicy,
	    AdjFwdSimPolicy,
	    AdjSimPolicy
	    > const & x)
      : dom(x.dom), rng(x.rng), 
	stream(x.stream), 
	minit(x.minit),
	pars(NULL),
	lpars(NULL),
	dump_steps(x.dump_steps),
	dump_pars(x.dump_pars),
	dump_term(x.dump_term)
    {
      // cerr<<"literal copy, so pararray is copied"<<endl;
      int err=0;
      if (err=ps_copy(&pars,*(x.pars))) {
	RVLException e;
	e<<"Error: IWOP copy constructor from ps_copy, err="<<err<<"\n";
	throw e;
      }
      /* see what we've got */
      if (dump_pars) {
	fprintf(stream,"PARS IN IWOP COPY CONSTRUCTOR\n");
	ps_printall(*pars,stream);
      }
    }

    ~IWaveOp() { 
      ps_delete(&pars);
      if (lpars) ps_delete(&lpars);
      /* cerr<<"good bye!\n"; */ }
      
    const Space<ireal> & getDomain() const { return dom; }
    const Space<ireal> & getRange() const { return rng; }

    // added 23.06.10 to facilitate using source as variable
    // without admitting that it's part of domain
    PARARRAY & getPar() { return *pars; }
    PARARRAY const & getPar() const { return *pars; }

    ostream & write(ostream & str) const {
      str<<"IWOP object\n";
      return str;
    }
  };
}


#endif
