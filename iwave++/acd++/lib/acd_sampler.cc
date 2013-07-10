#include "acd_sampler.hh"

#define ACD_SAM_VERB

namespace ACD {

  void cd_linreadmedia(RDOM dom,
		       FILE * stream,
		       PARARRAY par,
		       int panelindex) {
    /*==================================================
      ================= DECLARATIONS ===================
      ==================================================*/

    /* other workspace */
    IPNT ran;        /* axis lens */
    IPNT rags;       /* axis starts */
    IPNT rage;       /* axis ends */

    /* error flag */
    int err=0;
 
    /*==================================================
      =============== END DECLARATIONS =================
      ==================================================*/  

    /* initializations */

    /* axis lens, starts for csq mod - allocated, not virtual */
    rd_size(&dom,D_CSQ,ran);
    rd_gse(&dom, D_CSQ,rags,rage);

    /* FILE I/O SECTION */
    string dcsqname="";
    if (!parse<string>(par,"dcsq",dcsqname)) {
      RVLException e;
      e<<"Error: acd++::cd_linread\n";
      e<<"  failed to read key dcsq\n";
      throw e;
    }
		  
    err=rsfread(dom._s[D_CSQ]._s0, rags, ran, dcsqname.c_str(), 1, stream, panelindex); 
    if (err) {
      RVLException e;
      e<<"Error: acd++::cd_linread from rsfread, err="<<err<<"\n";
      e<<"  failed to read from file = "<<dcsqname<<"\n";
      throw e;    
    }
		      
  }

  ACDSampler::ACDSampler(IWaveState & _state)
    : state(_state), 
      tt(_state), 
      tstart(ct),
      tfinal(cte), 
      aps(_state,trace.t.tg),
      trs(_state,trace),
      is_samplerinit(false) {

    int err=0;
    IPNT tmp;
    for (int i=0;i<RARR_MAX_NDIM;i++) tmp[i]=D_UC;
    if ((err=sampler_construct(&trace,
			      &(state.getPAR()),
			      tmp,
			      RPNT_1,
			      RPNT_0,
			      0,
			      "hdrfile",
			      "datafile",
			       state.getStream()))) {
      RVLException e;
      e<<"Error: ACDSampler constructor from sampler_construct\n";
      e<<"returned code "<<err<<"\n";
      throw e;
    }
  }

  ACDSampler::~ACDSampler() {

    int err=0;

    if ((err=sampler_destroy(&trace))) {
      RVLException e;
      e<<"Error: ACDSampler destructor from sampler_destroy\n";
      e<<"returned code "<<err<<"\n";
      throw e;
    }
  }

  int ACDSampler::getNumRecs() const { return trace.t.tg.nrec; }
  int ACDSampler::getRecIndx() const {
    int panelindex;
    if (is_samplerinit) panelindex = trace.t.tg.irec;
    else panelindex = trace.t.tg.xrec;
    return panelindex;
  }  
  int ACDSampler::getFirstRec() const { return trace.t.tg.first; }
  int ACDSampler::getLastRec() const { return trace.t.tg.last; }

  bool ACDSampler::init() {
 
    int err = !(trace.t.tg.xrec < trace.t.tg.last+1);

    IMODEL & mdl = (state.IWaveState::getIWAVE()).model;
    FD_MODEL * fdm = (FD_MODEL *)(mdl.specs);

    if (!err) { 
      err=sampler_init(&trace,
		       &mdl,
		       &(state.getPAR()),
		       state.getStream());
    }
    // this should be the only terminal output in the entire run
#ifdef ACD_SAM_VERB
    if (!err) { 
      if ((retrieveGlobalRank()==0) && (trace.t.tg.irec == trace.t.tg.first)) 
	cerr<<"IWAVE++ Acoustic Constant Density\n";
      cerr<<	" fwd ->"<<
		" rkw="<<retrieveGlobalRank()<<
		" rkl="<<retrieveRank()<<
		" src="<<trace.t.tg.irec<<endl;
    }
#endif
    if (!err) is_samplerinit = true;

    // fprintf(state.getStream(),"ACDSampler::init, finish sampler_init, err =%d\n",err);
    if (!err) 
      aps.init();
    
    // fprintf(state.getStream(),"ACDSampler::init, finish aps.init(), err =%d\n",err);
    // cerr<< "ACDSampler::init, xrec = "<<trace.t.tg.xrec<<", nrec = "<<trace.t.tg.nrec<<", err = "<<err<<endl; 
    

    /* set start time */  
    if (!err) {
      TIMESTEPINDEX & t = tstart.getCstruct();
      t.it=aps.getStartIndex();
      //      cerr<<"start index = "<<t.it<<endl;
      t.iv=0;
      t.dt=trace.t.tg.dt;
    }

    /* set end time */
    if (!err) {
      TIMESTEPINDEX & te = tfinal.getCstruct();
      // te.it= trace.t.istop + 1;
      // te.iv= 0;
      te.it= trace.t.istop;
      te.iv = fdm->numsubsteps() - 1;
      te.dt= trace.t.tg.dt;

      TSIndex _tend(te);
      tt.setTargetTime(_tend);
      //     tt.setTargetTime(trace.t.istop + 1,trace.t.tg.dt);
      // previous version:  tt.setTargetTime(trace.t.tg.nt,trace.t.tg.dt);  
    }
    // cerr<< "exit ACDSampler::init(), with err = "<<err<<"\n";
    //cin>>err;
    if (err) return false;
    return true;

  }

  void ACDSampler::flush() const {

    /* extract step, origin vectors */
    RPNT d;
    RPNT o; 
    get_d(d, state.IWaveState::getIWAVE().model.g);
    get_o(o, state.IWaveState::getIWAVE().model.g);

    //    cerr<<"ACDSampler::flush: nt="<<trace.t.tg.nt<<" ntout="<<trace.t.tg.ntout<<" t0="<<trace.t.tg.t0<<" t0out="<<trace.t.tg.t0out<<endl;
    //    cerr<<"ACDSampler::flush: ptr to tg = "<<&(trace.t.tg)<<endl;

    /* write out traes to file */
    int err=writetraces(&(trace.t.tg),d,o,state.getStream());
    
    fflush(state.getStream());

    /* reset is_samplerinit to false for the next record */
    is_samplerinit = false;

    /*l
      fprintf(stderr,"at end of ACDSampler::flush\n");
      scanf("%d",&err);
    */
    if (err) {
      RVLException e;
      e<<"ERROR: ACDSampler::flush from tracegeom::writetraces. ABORT\n";
      throw e;
    }
    
  }

  ostream & ACDSampler::write(ostream & str) const {
    sampler_fprint(&trace,state.getStream());
    aps.fprint();
    return str;
  }
  
  /* ACDLinSampler */

  int ACDLinSampler::getNumRecs() const { return trace.t.tg.nrec; }
  int ACDLinSampler::getRecIndx() const {
    int panelindex;
    if (is_samplerinit) panelindex = trace.t.tg.irec;
    else panelindex = trace.t.tg.xrec;
    return panelindex;
  }  
  int ACDLinSampler::getFirstRec() const { return trace.t.tg.first; }
  int ACDLinSampler::getLastRec() const { return trace.t.tg.last; }

  ACDLinSampler::ACDLinSampler(IWaveLinState & _state)
    : state(_state), 
      tt(_state), 
      tstart(ct), 
      tfinal(cte),
      ltt(_state),
      aps(_state,trace.t.tg),
      trs(_state,trace),
      ltrs(_state,trace),
      //      abs(_state),
      is_samplerinit(false) {

    /***************************/
    //    cerr<<"beg ACDLinSampler\n";
    //    iwave_fprintall(stderr);
    /***************************/

    int err=0;    

    IPNT tmp;
    for (int i=0;i<RARR_MAX_NDIM;i++) tmp[i]=D_UC;			      
    if ((err=sampler_construct(&trace,
			      &(state.getPAR()),
			      tmp,
			      RPNT_1,
			      RPNT_0,
			      0,
			      "hdrfile",
			      "ldatafile",
			       state.getStream()))) {
      RVLException e;
      e<<"Error: TraceSampler constructor from sampler_construct\n";
      e<<"returned code "<<err<<"\n";
      throw e;
    }
    //    cerr<<"end ACDLinSampler\n";
  }

  ACDLinSampler::~ACDLinSampler() {
    int err=0;
    if ((err=sampler_destroy(&trace))) {
      RVLException e;
      e<<"Error: ACDSampler destructor from sampler_destroy\n";
      e<<"returned code "<<err<<"\n";
      throw e;
    }
  }
  
  bool ACDLinSampler::init() {
    try {

      int printact = state.getIWAVE().printact;
      FILE *stream = state.getStream();
      if(printact>1){
	fprintf(stream,"\n---enter ACDLinSampler::init()\n ");
      }
      int err = !(trace.t.tg.xrec < trace.t.tg.last+1);
      //    cerr<< "ACDLinSampler::init, xrec = "<<trace.t.tg.xrec<<", nrec = "<<trace.t.tg.nrec<<", trace.t.tg.last = "<<trace.t.tg.last<<", panelindex = "<<panelindex<<", err = "<<err<<endl;

      // iwave_lininit contains: 
      //    (1) call iwave_init to set up the reference field, i.e.:
      //      (a) zero dynamic fields
      //      (b) call mread func (input ref model, update time step dt and tspars) 
      //    (2) set up the perturbed field (another IWAVE object), i.e.:
      //      (a) zero dynamic fields   
      //      (b) make sure that the pointers of non-dynamic fields point to ref model
      //      (c) copy dt and tspars from ref model   
      //      (d) call mlinread to read perterbed coeff into the non-dynamic fields                       
    
      /* linearized simulation needs same "specs" as reference
	 same with tspars and post-time-step (boundary condns) */
      /* HYPOTHESIS: this is all redundant now - use pts from ref state */

      IMODEL * refmodel = &(state.IWaveState::getIWAVE().model);
      IMODEL * linmodel = &(state.getIWAVE().model);

      FD_MODEL * rfdm = (FD_MODEL *) (refmodel->specs);
      FD_MODEL * lfdm = (FD_MODEL *) (linmodel->specs);
    
      // if linearized model is extended (whether ref model is, or not)
      // or if this is the first record sim in the block, load linearized 
      // model
      if ( (linmodel->g.dim < linmodel->g.gdim) || 
	   (trace.t.tg.xrec == trace.t.tg.first)) {	
	/* start of static_init block */
	
	if (!err){
	
	  ACD_TS_PARS * reftspars = (ACD_TS_PARS *)(rfdm->fdpars);
	  ACD_TS_PARS * lintspars = (ACD_TS_PARS *)(lfdm->fdpars);
	  // copy dt and tspars from ref model to pert model
	  rfdm->parcopy(lintspars,reftspars);
	  (linmodel->tsind).dt = (refmodel->tsind).dt;
	}
      
	/*-----------------------------------------------------------------*  
	 * set perturbed buoyancy and dcsq modulus -----------------------*
	 *-----------------------------------------------------------------*/
	//cerr<<"ACDLinSampler::init, calling cd_linreadmedia, err ="<<err<<endl;
	//cin >> err;

	if (!err) 
	  cd_linreadmedia(linmodel->ld_a,
			  state.getStream(),
			  state.getPAR(),
			  this->getRecIndx());
      } /* end of static_init block */
    
      //cerr<<"ACDLinSampler::init, calling sampler_init, err ="<<err<<endl;
    
      if (!err)
	err=sampler_init(&trace,
			 linmodel,
			 &(state.getPAR()),
			 state.getStream());
    
      // this should be the only terminal output in the entire run
#ifdef ACD_SAM_VERB
      if (!err) {
	if ((retrieveGlobalRank()==0) && (trace.t.tg.irec == trace.t.tg.first)) 
	  cerr<<"IWAVE++ Acoustic Constant Density\n";
	cerr<<	  "lin ->"<<
		  " rkw="<<retrieveGlobalRank()<<
		  " rkl="<<retrieveRank()<<
		  " src="<<trace.t.tg.irec<<endl;
      }
#endif
      //cerr<<"ACDLinSampler::init, calling aps.init, err ="<<err<<endl;
      //cin >> err;
      if (!err)    aps.init();
   
      // set start time   
      if (!err) {
	TIMESTEPINDEX & t = tstart.getCstruct();
	t.it=aps.getStartIndex();
	t.iv=0;
	t.dt=trace.t.tg.dt;
      }

      // set end time 
      if (!err) { 
	TIMESTEPINDEX & te = tfinal.getCstruct();
	//te.it= trace.t.istop + 1;
	// te.iv=0;
	te.it= trace.t.istop;
	te.iv = rfdm->numsubsteps() - 1;
	te.dt= trace.t.tg.dt;

	TSIndex _tgt(te);
	//      fprintf(stderr,"in ACDLinSampler::init: nt=%d istop=%d dt=%e\n",trace.t.tg.nt,trace.t.istop,trace.t.tg.dt);
	// tt.setTargetTime(trace.t.istop + 1,trace.t.tg.dt);
	tt.setTargetTime(_tgt);
	// previous version:  tt.setTargetTime(trace.t.tg.nt,trace.t.tg.dt); 
	//  ltt.setTargetTime(trace.t.istop + 1,trace.t.tg.dt); 
	ltt.setTargetTime(_tgt);      
      }
  
      if(printact>1){
	fprintf(stream,"\n---exit ACDLinSampler::init()\n ");
      }
    
      //cerr<<"ACDLinSampler::init, exiting, err ="<<err<<endl;
      //cin >> err;
      if (err) return false;
      return true;
    }
    catch (RVLException & e) {
      e<<"\ncalled from ACDLinSampler::init\n";
      throw e;
    }
  }

  void ACDLinSampler::flush() const {

    // extract step, origin vectors 
    RPNT d;
    RPNT o; 
    get_d(d, state.getIWAVE().model.g);
    get_o(o, state.getIWAVE().model.g);

    // write out traes to file 
    int err=writetraces(&(trace.t.tg),d,o,state.getStream());
    if (err) {
      RVLException e;
      e<<"ERROR: ACDLinSampler::flush from tracegeom::writetraces. ABORT \n";
      throw e;
    }
    //cerr<<"after writetraces \n";

    /* reset is_samplerinit to false for the next record */
    is_samplerinit = false;

  } 

  ostream & ACDLinSampler::write(ostream & str) const {
    sampler_fprint(&trace,state.getStream());
    aps.fprint();
    return str;
  }
 
  /* ACDAdjSampler */
 
  ACDAdjSampler::ACDAdjSampler(IWaveLinState & _state)
    : state(_state), 
      tt(_state,false,cerr), 
      bwdtt(_state,false,cerr), 
      tstart(ct), 
      tfinal(cte),
      aps(_state,trace.t.tg),
      trs(_state,trace),
      adjtrs(_state,trace),
      is_samplerinit(false) {  

    int err=0;

    IPNT tmp;
    for (int i=0;i<RARR_MAX_NDIM;i++) tmp[i]=D_UC;

    // mod 02.13: initbuf=-1 causes adjoint spline interp, as opposed to
    // initbuf=1 which is straight spline interp, onto internal time grid
    if ((err=sampler_construct(&trace,
			       &(state.getPAR()),
			       tmp,
			       RPNT_1,
			       RPNT_0,
			       1,
			       "hdrfile",
			      "rdatafile",
			       state.getStream()))) {
      RVLException e;
      e<<"Error: ACDAdjSampler constructor from sampler_construct\n";
      e<<"returned code "<<err<<"\n";
      throw e;
    }

    // adjoint sampler requires that number of panels (separate output grids)
    // be equal to number of records (for data-domain extended modeling) or to 
    // 1 (for non-extended modeling)
    if ((get_panelnum_grid(state.getIWAVE().model.g) != 1) && 
	(get_panelnum_grid(state.getIWAVE().model.g) != trace.t.tg.nrec)) {
      RVLException e;
      e<<"Error: ACDAdjSampler constructor\n";
      e<<"number of model output bins must equal 1 (standard modeling) \n";
      e<<"or number of data records = "<<trace.t.tg.nrec<<", however = ";
      e<<get_panelnum_grid(state.getIWAVE().model.g) <<"\n";
      throw e;
    }

  }

  ACDAdjSampler::~ACDAdjSampler() {
   
    int err=0;

    if ((err=sampler_destroy(&trace))) {
      RVLException e;
      e<<"Error: ACDSampler destructor from sampler_destroy\n";
      e<<"returned code "<<err<<"\n";
      throw e;
    }
  }
  
  int ACDAdjSampler::getNumRecs() const { return trace.t.tg.nrec; }
  int ACDAdjSampler::getRecIndx() const {
    int panelindex;
    if (is_samplerinit) panelindex = trace.t.tg.irec;
    else panelindex = trace.t.tg.xrec;
    return panelindex;
  }  
  int ACDAdjSampler::getFirstRec() const { return trace.t.tg.first; }
  int ACDAdjSampler::getLastRec() const { return trace.t.tg.last; }

  bool ACDAdjSampler::init() {

    int printact = state.getIWAVE().printact;
    FILE *stream = state.getStream();
    int err = !(trace.t.tg.xrec < trace.t.tg.last+1);
    
    if (printact) fprintf(stream,"ACDAdjSampler::init, xrec=%d, nrec=%d last=%d err=%d\n",trace.t.tg.xrec,trace.t.tg.nrec,trace.t.tg.last,err);

    trace.t.index = D_UC;  

    IMODEL * fwdmodel = &(state.IWaveState::getIWAVE().model);
    IMODEL * bwdmodel = &(state.getIWAVE().model);
    
    FD_MODEL * fwdfdm = (FD_MODEL *) (fwdmodel->specs);
    FD_MODEL * bwdfdm = (FD_MODEL *) (bwdmodel->specs);

    if ((trace.t.tg.xrec == trace.t.tg.first) ||
	(get_panelnum_grid(state.getIWAVE().model.g) > 1)) {

      // for first record in group, if output record has not been
      // visited before, zero out output rarrays (D_P and D_MV of 
      // pert domain)
      // Expls:
      // if get_panelnum_grid = 1 (standard modeling), then quot.quot = 
      // panelindex, and is only zero irec=0. Therefore output is zeroed
      // only for the first record of the first group, and accumulates for
      // all others
      // if get_panelnum_grid = nrec, then quot.quot = 0 in all cases,
      // as input and output records are in 1-1 correspondence, so output
      // should be zeroed for every record.
      // Note that output is reduced over the group for the first case
      /* start of static_init block */
      if (!err) {
	//backward simulation needs same "specs" as reference
	// same with tspars and post-time-step (boundary condns) 
	
	ACD_TS_PARS * fwdtspars = (ACD_TS_PARS *)(fwdfdm->fdpars);
	ACD_TS_PARS * bwdtspars = (ACD_TS_PARS *)(bwdfdm->fdpars);
	// copy dt and tspars from ref model to pert model
	fwdfdm->parcopy(bwdtspars,fwdtspars);
	(bwdmodel->tsind).dt = (fwdmodel->tsind).dt;
	
	// initializing  mig_model 
	RDOM * ddom = &(bwdmodel->ld_a);
	ra_a_zero(&(ddom->_s[D_CSQ]));
	
      }// end of 'if(!err)...' branch 
    } /* end of static_init block */
    
    // fprintf(stream,"ACDAdjSampler::init, before sampler_init, err =%d\n",err);
    //    cerr<<"ACDAdjSampler::init, before sampler_init, err = "<<err<<endl;
    if (!err)
      err=sampler_init(&trace,
		       bwdmodel,
		       &(state.getPAR()),
		       state.getStream());

    // this should be the only terminal output in the entire run
#ifdef ACD_SAM_VERB
    if (!err) {
      if ((retrieveGlobalRank()==0) && (trace.t.tg.irec == trace.t.tg.first)) 
	cerr<<"IWAVE++ Acoustic Constant Density\n";
      cerr<<    "adj ->"<<
		" rkw="<<retrieveGlobalRank()<<
		" rkl="<<retrieveRank()<<
		" src="<<trace.t.tg.irec<<endl;
    }
#endif
    if(!err) is_samplerinit = true;
    
    //cerr<<"ACDAdjSampler::init, after sampler_init, err = "<<err<<endl;

    //fprintf(stream,"ACDAdjSampler::init, before aps.init, err =%d\n",err);
    fflush(stream);
    //cerr<<"ACDAdjSampler::init, before aps.init, err = "<<err<<endl;
    if (!err)
      aps.init();
   
    //cerr<<"ACDAdjSampler::init, after aps.init, err = "<<err<<endl;
    
    /* set start and final time*/   
    if (!err) {
      TIMESTEPINDEX & t = tstart.getCstruct();
      t.it=aps.getStartIndex();
      t.iv=0;
      t.dt=trace.t.tg.dt;
    
      TIMESTEPINDEX & te = tfinal.getCstruct();
      te.it= trace.t.istop;
      te.iv = fwdfdm->numsubsteps() - 1;
      te.dt= trace.t.tg.dt;
      
      // set target time for  forward and backward terminators 
      // fprintf(stderr,"in ACDAdjSampler::init: nt=%d istop=%d dt=%e\n",trace.t.tg.nt,trace.t.istop,trace.t.tg.dt);
      //  tt.setTargetTime(trace.t.istop,trace.t.tg.dt);
      //  bwdtt.setTargetTime(trace.t.istart,trace.t.tg.dt);
      TSIndex _tstart(t);
      TSIndex _tend(te);
      tt.setTargetTime(_tend);
      bwdtt.setTargetTime(_tstart);
    }
  
    //    fprintf(stream,"\n---exit ACDAdjSampler::init()\n ");
    //    fflush(stream);
   
    if (printact) fprintf(stream,"ACDAdjSampler::init, exit with err = %d\n",err);
     
    if (err) return false;
    return true;
  }

  void ACDAdjSampler::flush() const {

    // this is a no-op if model is physical (not extended)
    // and simulation is not finished
    // if physical model, invoke barrier to ensure that all contributions 
    // are available to reduction
    FILE *stream = state.getStream(); 
    fflush(stream);
    if (get_panelnum_grid(state.getIWAVE().model.g) == 1) {
      if (this->getRecIndx() < this->getLastRec()) return;
#ifdef IWAVE_USE_MPI
      MPI_Barrier(retrieveGlobalComm());
#endif
    }
    PARARRAY & pars = state.getPAR();

    // fish output filename out of param table
    string mcsqname="";
    if (!parse<string>(pars,"mcsq",mcsqname)) {
      RVLException e;
      e<<"Error: ACDAdjSampler::flush\n";
      e<<"  failed to read mcsqmod from par file\n";
      throw e;
    }

    int err = 0;
    int dim;
    IPNT rags, rage;
    IPNT ran;
        
    rd_ndim(&(state.getIWAVE().model.ld_a),D_CSQ,&dim);
    rd_size(&(state.getIWAVE().model.ld_a),D_CSQ,ran);
    rd_gse(&(state.getIWAVE().model.ld_a),D_CSQ,rags,rage);

    // for non-extended models, reduce group adjoint accumulations to global
#ifdef IWAVE_USE_MPI
    /*
    cerr<<"----------ACDAdjSampler::flush csqmod="<<migcsqname
	<<" rkw="<<retrieveGlobalRank()
	<<" rkl="<<retrieveRank()
	<<" gid="<<retrieveGroupID()
	<<" pnm="<<get_panelnum_grid(state.getIWAVE().model.g)
	<<endl;
    */
    if (get_panelnum_grid(state.getIWAVE().model.g) == 1) {
      // allocate buffer for reduction target
      int ntot=1;
      for (int i=0;i<dim;i++) ntot*=ran[i];
      ireal * buf = new ireal[ntot];
      // reduce along remote comm axis - same domain in every group
      err=MPI_Reduce(state.getIWAVE().model.ld_a._s[D_CSQ]._s0,buf,ntot,IWAVE_MPI_REAL,MPI_SUM,0,retrieveRemComm());
      if (err) {
	RVLException e;
	e<<"Error: AcdAdjSampler::flush() from MPI_Reduce, D_CSQ\n";
	throw e;
      }
      // on group 0, copy back to alloca domain
      if (retrieveGroupID()==0) memcpy(state.getIWAVE().model.ld_a._s[D_CSQ]._s0,buf,ntot*sizeof(ireal));
      delete [] buf;
    }
   
    int g1 = get_panelnum_grid(state.getIWAVE().model.g);
    int g2 = retrieveGroupID();
    // write occurs only if (i) extended modeling - all groups, or (ii) group 0, otherwise
    if ((g1 > 1) || (g2==0)) {

#endif

      //      cerr<<"---------ACDAdjSampler::flush - write csqmod gid="<<retrieveGroupID()<<endl; 
      //      fprintf(state.getStream(),"---------ACDAdjSampler::flush migrated csqmod = %s record=%d\n",migcsqname,this->getRecIndx());
      err = rsfwrite(state.getIWAVE().model.ld_a._s[D_CSQ]._s0,
		     rags,
		     ran,
		     mcsqname.c_str(),
		     0,
		     state.getStream(),
		     this->getRecIndx()); 
    
      //      cerr<<"-----------ACDAdjSampler::flush, wrote migrated csq-modulus to "<<migcsqname<<", err = "<<err<<endl;
      if (err) {
	RVLException e;
	e<<"Error: ACDAdjSampler::flush()\n";
	e<<"rsfwrite return err for writing mig-dkappa\n";
	throw e;
      }

#ifdef IWAVE_USE_MPI
    }
#endif

    /* reset is_samplerinit to false for next record */
    is_samplerinit = false;

  } 


  ostream & ACDAdjSampler::write(ostream & str) const {
    sampler_fprint(&trace,state.getStream());
    aps.fprint();
    return str;
  }

}

