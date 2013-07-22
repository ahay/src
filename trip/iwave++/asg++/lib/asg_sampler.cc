#include "asg_sampler.hh"

/** Input perturbed parameters.
    Target: perturbed bulk modulus and buoyancy on appropriately shifted grids.
    Set REAL_ZERO as default values. 

    D.S. 09.06.09  
*/

int sg_linreadmedia(RDOM dom,
		    FILE * stream,
		    PARARRAY par,
		    int panelindex){
  /*==================================================
    ================= DECLARATIONS ===================
    ==================================================*/

  /* workspace for filenames */  
  char * dkappakey=NULL;
  char * dbulkkey=NULL;
  char * dbuoykey=NULL;

  /* default values */
  ireal refdkappa = REAL_ZERO;
  ireal refdbuoy = REAL_ZERO;

  /* tmp storage for read */
  RARR Btmp;

  /* other workspace */
  int i,j,m; /* counters */
  int dim;         /* problem dimension */
  IPNT ran;        /* axis lens */
  IPNT rags;       /* axis starts */
  IPNT rage;       /* axis ends */
  unsigned long ntot;        /* total num of words in phys dom */

  /* error flag */
  int err=0;
  IPNT ip;
  IPNT raVgs, raVge;
  ireal q;
  int mult;
 
  /*==================================================
    =============== END DECLARATIONS =================
    ==================================================*/  

  /* initializations */

  ra_setnull(&Btmp);
  /* axis lens, starts for bulk mod - allocated, not virtual */
  rd_ndim(&dom,D_MP0,&dim);
  rd_size(&dom,D_MP0,ran);
  rd_gse(&dom, D_MP0,rags,rage);
  ntot=1;
  for (i=0;i<dim;i++) ntot*=ran[i];

  /* SANITY CHECK MOVED TO SINGLE BOYANCY - I&T 03/05 */

  /* FILE I/O SECTION */

  /* read dbulkmod into D_MP0 */
  if (!ps_flcstring(par,"dkappa", &dkappakey) ) {
    if (!ps_flcstring(par,"dbulkmod", &dbulkkey)) {
      userfree_(dkappakey);
      userfree_(dbulkkey);
      fprintf(stream,
	      "Error: sg_linreadmedia from ra_rsfread: both dkappa and dbulkmod"
	      " are supplied (ambiguity)\n");
      return E_BADINPUT;
    }
		  
    err=rsfread(dom._s[D_MP0]._s0, rags, ran, dkappakey, 1, stream, panelindex); 
    
    userfree_(dkappakey);
    if (err) {
      fprintf(stream,"Error: sg_linreadmedia from ra_rsfread - dkappa\n");
      return err;
    }
  }
  else if (!ps_flcstring(par,"dbulkmod", &dbulkkey)) {
    err=rsfread(dom._s[D_MP0]._s0, rags, ran, dbulkkey, 1, stream, panelindex);
    userfree_(dbulkkey);
    if (err) {
      fprintf(stream,"Error: sg_readmedia from ra_rsfread - dbulkmod\n");
      return err;
    }
  }
  else if (!ps_flcstring(par,"dvelocity", &dbulkkey)) {
    err=rsfread(dom._s[D_MP0]._s0, rags, ran, dbulkkey, 1, stream, panelindex);
    userfree_(dbulkkey);
    if (err) {
      fprintf(stream,"Error: sg_readmedia from ra_rsfread - dvelocity mod\n");
      return err;
    }
  }
  else {
    /* if no data read, set to ref value */
      for (j=0;j<(int)ntot;j++) dom._s[D_MP0]._s0[j]=refdkappa;
  }

  /* read dbuoyancy */
  i=0; 
  mult=0;
  /* read single field into workspace, create shifted arrays by averaging. */

  if (!mult) {
    rd_gse(&dom, D_MP0, rags, rage);
    /* extend bulk borders if density is outside */
    for (i=0; i<dim; i++) {
      rd_gse(&dom,D_MV[i],raVgs,raVge);
      for (j=0; j<dim; j++) {
	if (j == i) {
	  if ( raVgs[j] < rags[j] ) rags[j] = raVgs[j];
	  if ( raVge[j] >=rage[j] ) rage[j] = raVge[j] + 1;
	}
	else {
	  if ( raVgs[j] < rags[j] ) rags[j] = raVgs[j];
	  if ( raVge[j] > rage[j] ) rage[j] = raVge[j];
	}			
      }  
    }
    err = ra_create(&Btmp, dim, rags, rage);
    if ( err ) {
      fprintf(stream,
	      "Error: sg_linreadmedia from rsfread - cannot allocate tmp array.\n");
      return err;
    }
    ra_size(&Btmp, ran);
	  
    dbuoykey=NULL;


    if (!ps_flcstring(par,"dbuoyancy",&dbuoykey)) {
      err=rsfread(Btmp._s, rags, ran, dbuoykey, 1, stream, panelindex); 
      if (err) {
	fprintf(stream,
		"Error: sg_linreadmedia from rsfread - dbuoykey (single) = %s\n",dbuoykey);
	ra_destroy(&Btmp);
	return err;
      }

      /* shift along coord axes, average */
      for (i = 0; i < dim; i++) {
	rd_gse(&dom, D_MV[i], raVgs, raVge);
	for ( j = dim; j < RARR_MAX_NDIM; ++j ) 
	  raVgs[j] = raVge[j] = 0; /* extra dimensions */

#if RARR_MAX_NDIM > 2
	for ( ip[2] = raVgs[2]; ip[2] <= raVge[2]; ++ip[2] )
#endif
#if RARR_MAX_NDIM > 1
	  for ( ip[1] = raVgs[1]; ip[1] <= raVge[1]; ++ip[1] )
#endif
	    for ( ip[0] = raVgs[0]; ip[0] <= raVge[0]; ++ip[0] ) {
	      q = ra_gget(&Btmp, ip);
	      ip[i] += 1;
	      q += ra_gget(&Btmp, ip);
	      ip[i] -= 1;
	      rd_gset(&dom, D_MV[i], ip, q * 0.5);  /* cause segmentation fault */ 
	    }
      }	  
    }
    else if (!ps_flcstring(par,"ddensity",&dbuoykey)) {
      err=rsfread(Btmp._s, rags, ran, dbuoykey, 1, stream, panelindex); 
      if (err) {
	fprintf(stream,
		"Error: sg_linreadmedia from rsfread - (ddensity) dbuoykey (single) = %s\n",dbuoykey);
	ra_destroy(&Btmp);
	return err;
      }

      /* shift along coord axes, average */
      for (i = 0; i < dim; i++) {
	rd_gse(&dom, D_MV[i], raVgs, raVge);
	for ( j = dim; j < RARR_MAX_NDIM; ++j ) 
	  raVgs[j] = raVge[j] = 0; /* extra dimensions */

#if RARR_MAX_NDIM > 2
	for ( ip[2] = raVgs[2]; ip[2] <= raVge[2]; ++ip[2] )
#endif
#if RARR_MAX_NDIM > 1
	  for ( ip[1] = raVgs[1]; ip[1] <= raVge[1]; ++ip[1] )
#endif
	    for ( ip[0] = raVgs[0]; ip[0] <= raVge[0]; ++ip[0] ) {
	      q = ra_gget(&Btmp, ip);
	      ip[i] += 1;
	      q += ra_gget(&Btmp, ip);
	      ip[i] -= 1;
	      rd_gset(&dom, D_MV[i], ip, q * 0.5);  /* cause segmentation fault */ 
	    }
      }	  
    }
    else {
      /* last resort: default value */
      for (i=0;i<dim;i++) {
	rd_size(&dom,D_MV[i],ran);
	m=1;
	for (j=0;j<dim;j++) m*=ran[j];
	for (j=0;j<m;j++) dom._s[D_MV[i]]._s0[j]=refdbuoy;
      }
    }
   
    ra_destroy(&Btmp);
  
  }

  return err;   
		      
}

namespace ASG {

  ASGSampler::ASGSampler(IWaveState & _state)
    : state(_state), 
      trs(_state,trace),
      aps(_state,trace.t.tg),
      tt(_state), 
      tstart(ct),
      tfinal(cte), 
      is_samplerinit(false) {
   
    //    ps_printall(state.getPAR(),stderr);
    int err=0;    
    RPNT mult;
    for (int i=0;i<RARR_MAX_NDIM;i++) mult[i]=REAL_ONE/((ireal)(state.getIWAVE().model.g.dim));

    if ((err=sampler_construct(&trace,
			      &(state.getPAR()),
			      D_P,
			      mult,
			      RPNT_0,
			      0,
			      "hdrfile",
			      "datafile",
			       state.getStream()))) {
      RVLException e;
      e<<"Error: TraceSampler constructor from sampler_construct\n";
      e<<"returned code "<<err<<"\n";
      throw e;
    }
    else {
      fprintf(state.getStream(),"NOTE: TraceSampler constructed\n");
    }
  }

  ASGSampler::~ASGSampler() {

    int err=0;

    if ((err=sampler_destroy(&trace))) {
      RVLException e;
      e<<"Error: ASGSampler destructor from sampler_destroy\n";
      e<<"returned code "<<err<<"\n";
      throw e;
    }
  }

  int ASGSampler::getNumRecs() const { return trace.t.tg.nrec; }
  int ASGSampler::getRecIndx() const {
    int panelindex;
    if (is_samplerinit) panelindex = trace.t.tg.irec;
    else panelindex = trace.t.tg.xrec;
    return panelindex;
  }  
  int ASGSampler::getFirstRec() const { return trace.t.tg.first; }
  int ASGSampler::getLastRec() const { return trace.t.tg.last; }

  bool ASGSampler::init() {
 
    // int err=!(trace.t.tg.xrec < trace.t.tg.nrec);
    // fprintf(state.getStream(),"ASGSampler::init, xrec = %d nrec=%d err=%d\n",trace.t.tg.xrec,trace.t.tg.nrec,err);
    int err = !(trace.t.tg.xrec < trace.t.tg.last+1);
    //    cerr<< "ASGSampler::init, xrec = "<<trace.t.tg.xrec<<", nrec = "<<trace.t.tg.nrec<<", trace.t.tg.last = "<<trace.t.tg.last<<", err = "<<err<<endl;   

    IMODEL & mdl = (state.IWaveState::getIWAVE()).model;
    FD_MODEL * fdm = (FD_MODEL *)(mdl.specs);

    if (!err) 

      err=sampler_init(&trace,
		       &mdl,
		       &(state.getPAR()),
		       state.getStream());

    // this should be the only terminal output in the entire run
    if (!err) cerr<<
                "fwd ->"<<
		" rkw="<<retrieveGlobalRank()<<
		" rkl="<<retrieveRank()<<
		" src="<<trace.t.tg.irec<<endl;

    if (!err) is_samplerinit = true;

    // fprintf(state.getStream(),"ASGSampler::init, finish sampler_init, err =%d\n",err);
    if (!err) 
      aps.init();
    
    // fprintf(state.getStream(),"ASGSampler::init, finish aps.init(), err =%d\n",err);
    // cerr<< "ASGSampler::init, xrec = "<<trace.t.tg.xrec<<", nrec = "<<trace.t.tg.nrec<<", err = "<<err<<endl; 
    

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
    // cerr<< "exit ASGSampler::init(), with err = "<<err<<"\n";
    //cin>>err;
    if (err) return false;
    return true;

  }

  void ASGSampler::flush() const {

    /* extract step, origin vectors */
    RPNT d;
    RPNT o; 
    get_d(d, state.IWaveState::getIWAVE().model.g);
    get_o(o, state.IWaveState::getIWAVE().model.g);

    //    cerr<<"ASGSampler::flush: nt="<<trace.t.tg.nt<<" ntout="<<trace.t.tg.ntout<<" t0="<<trace.t.tg.t0<<" t0out="<<trace.t.tg.t0out<<endl;
    //    cerr<<"ASGSampler::flush: ptr to tg = "<<&(trace.t.tg)<<endl;

    /* write out traes to file */
    int err=writetraces(&(trace.t.tg),d,o,state.getStream());
    
    fflush(state.getStream());

    /* reset is_samplerinit to false for the next record */
    is_samplerinit = false;

    /*l
      fprintf(stderr,"at end of ASGSampler::flush\n");
      scanf("%d",&err);
    */
    if (err) {
      RVLException e;
      e<<"ERROR: ASGSampler::flush from tracegeom::writetraces. ABORT\n";
      throw e;
    }
    
  }

  ostream & ASGSampler::write(ostream & str) const {
    sampler_fprint(&trace,state.getStream());
    aps.fprint();
    return str;
  }
  
  /* ASGLinSampler */

  int ASGLinSampler::getNumRecs() const { return trace.t.tg.nrec; }
  int ASGLinSampler::getRecIndx() const {
    int panelindex;
    if (is_samplerinit) panelindex = trace.t.tg.irec;
    else panelindex = trace.t.tg.xrec;
    return panelindex;
  }  
  int ASGLinSampler::getFirstRec() const { return trace.t.tg.first; }
  int ASGLinSampler::getLastRec() const { return trace.t.tg.last; }

  ASGLinSampler::ASGLinSampler(IWaveLinState & _state)
    : state(_state), 
      trs(_state,trace),
      ltrs(_state,trace),
      aps(_state,trace.t.tg),
      tt(_state), 
      ltt(_state),
      tstart(ct), 
      tfinal(cte),
      //      abs(_state),
      is_samplerinit(false) {
     
    int err=0;    

    RPNT mult;
    for (int i=0;i<RARR_MAX_NDIM;i++) mult[i]=REAL_ONE/((ireal)(state.getIWAVE().model.g.dim));

    if ((err=sampler_construct(&trace,
			      &(state.getPAR()),
			      D_P,
			      mult,
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

  }

  ASGLinSampler::~ASGLinSampler() {
    int err=0;
    if ((err=sampler_destroy(&trace))) {
      RVLException e;
      e<<"Error: ASGSampler destructor from sampler_destroy\n";
      e<<"returned code "<<err<<"\n";
      throw e;
    }
  }
  
  bool ASGLinSampler::init() {    
    int printact = state.getIWAVE().printact;
    FILE *stream = state.getStream();
    if(printact>1){
      fprintf(stream,"\n---enter ASGLinSampler::init()\n ");
    }
    int err = !(trace.t.tg.xrec < trace.t.tg.last+1);
    //    cerr<< "ASGLinSampler::init, xrec = "<<trace.t.tg.xrec<<", nrec = "<<trace.t.tg.nrec<<", trace.t.tg.last = "<<trace.t.tg.last<<", panelindex = "<<panelindex<<", err = "<<err<<endl;

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
	
	SGN_TS_PARS * reftspars = (SGN_TS_PARS *)(rfdm->fdpars);
	SGN_TS_PARS * lintspars = (SGN_TS_PARS *)(lfdm->fdpars);
	// copy dt and tspars from ref model to pert model
	rfdm->parcopy(lintspars,reftspars);
	(linmodel->tsind).dt = (refmodel->tsind).dt;
      }
      
      /*-----------------------------------------------------------------*  
       * set perturbed buoyancy and dbulk modulus -----------------------*
       *-----------------------------------------------------------------*/
      //cerr<<"ASGLinSampler::init, calling sg_linreadmedia, err ="<<err<<endl;
      //cin >> err;

      if (!err) 
	err = sg_linreadmedia(linmodel->ld_a,
			      state.getStream(),
			      state.getPAR(),
			      this->getRecIndx());
    } /* end of static_init block */
    
    //cerr<<"ASGLinSampler::init, calling sampler_init, err ="<<err<<endl;
    
    if (!err)
      err=sampler_init(&trace,
		       linmodel,
		       &(state.getPAR()),
		       state.getStream());
    
    // this should be the only terminal output in the entire run
    if (!err) cerr<<
                "lin ->"<<
		" rkw="<<retrieveGlobalRank()<<
		" rkl="<<retrieveRank()<<
		" src="<<trace.t.tg.irec<<endl;
    //cerr<<"ASGLinSampler::init, calling aps.init, err ="<<err<<endl;
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
      //      fprintf(stderr,"in ASGLinSampler::init: nt=%d istop=%d dt=%e\n",trace.t.tg.nt,trace.t.istop,trace.t.tg.dt);
      // tt.setTargetTime(trace.t.istop + 1,trace.t.tg.dt);
      tt.setTargetTime(_tgt);
      // previous version:  tt.setTargetTime(trace.t.tg.nt,trace.t.tg.dt); 
      //  ltt.setTargetTime(trace.t.istop + 1,trace.t.tg.dt); 
      ltt.setTargetTime(_tgt);      
    }
  
    if(printact>1){
      fprintf(stream,"\n---exit ASGLinSampler::init()\n ");
    }
    
    //cerr<<"ASGLinSampler::init, exiting, err ="<<err<<endl;
    //cin >> err;
    if (err) return false;
    return true;
  }

  void ASGLinSampler::flush() const {

    // extract step, origin vectors 
    RPNT d;
    RPNT o; 
    get_d(d, state.getIWAVE().model.g);
    get_o(o, state.getIWAVE().model.g);

    // write out traes to file 
    int err=writetraces(&(trace.t.tg),d,o,state.getStream());
    if (err) {
      RVLException e;
      e<<"ERROR: ASGLinSampler::flush from tracegeom::writetraces. ABORT \n";
      throw e;
    }
    //cerr<<"after writetraces \n";

    /* reset is_samplerinit to false for the next record */
    is_samplerinit = false;

  } 

  ostream & ASGLinSampler::write(ostream & str) const {
    sampler_fprint(&trace,state.getStream());
    aps.fprint();
    return str;
  }
 
  /* ASGAdjSampler */
 
  ASGAdjSampler::ASGAdjSampler(IWaveLinState & _state)
    : state(_state), 
      trs(_state,trace),
      adjtrs(_state,trace),
      aps(_state,trace.t.tg),
      tt(_state,false,cerr), 
      bwdtt(_state,false,cerr), 
      tstart(ct), 
      tfinal(cte),
      is_samplerinit(false) {  

    int err=0;
  
    /* note that this is a load op, i.e. an adjoint sampler */
    RPNT mult;
    for (int i=0;i<RARR_MAX_NDIM;i++) mult[i]=REAL_ONE/((ireal)(state.getIWAVE().model.g.dim));
    // mod 02.13: initbuf=-1 causes adjoint spline interp, as opposed to
    // initbuf=1 which is straight spline interp, onto internal time grid    
    if ((err=sampler_construct(&trace,
			      &(state.getPAR()),
			      D_P,
			      mult,
			      RPNT_0,
			      -1,
			      "hdrfile",
			      "rdatafile",
			       state.getStream()))) {
      RVLException e;
      e<<"Error: ASGAdjSampler constructor from sampler_construct\n";
      e<<"returned code "<<err<<"\n";
      throw e;
    }

    // adjoint sampler requires that number of panels (separate output grids)
    // be equal to number of records (for data-domain extended modeling) or to 
    // 1 (for non-extended modeling)
    if ((get_panelnum_grid(state.getIWAVE().model.g) != 1) && 
	(get_panelnum_grid(state.getIWAVE().model.g) != trace.t.tg.nrec)) {
      RVLException e;
      e<<"Error: ASGAdjSampler constructor\n";
      e<<"number of model output bins must equal 1 (standard modeling) \n";
      e<<"or number of data records = "<<trace.t.tg.nrec<<", however = ";
      e<<get_panelnum_grid(state.getIWAVE().model.g) <<"\n";
      throw e;
    }

  }

  ASGAdjSampler::~ASGAdjSampler() {
   
    int err=0;

    if ((err=sampler_destroy(&trace))) {
      RVLException e;
      e<<"Error: ASGSampler destructor from sampler_destroy\n";
      e<<"returned code "<<err<<"\n";
      throw e;
    }
  }
  
  int ASGAdjSampler::getNumRecs() const { return trace.t.tg.nrec; }
  int ASGAdjSampler::getRecIndx() const {
    int panelindex;
    if (is_samplerinit) panelindex = trace.t.tg.irec;
    else panelindex = trace.t.tg.xrec;
    return panelindex;
  }  
  int ASGAdjSampler::getFirstRec() const { return trace.t.tg.first; }
  int ASGAdjSampler::getLastRec() const { return trace.t.tg.last; }

  bool ASGAdjSampler::init() {

    int printact = state.getIWAVE().printact;
    FILE *stream = state.getStream();
    int err = !(trace.t.tg.xrec < trace.t.tg.last+1);
    
    if (printact) fprintf(stream,"ASGAdjSampler::init, xrec=%d, nrec=%d last=%d err=%d\n",trace.t.tg.xrec,trace.t.tg.nrec,trace.t.tg.last,err);

    trace.t.index = D_P[0];  // 10.10.10. (D.S. init for multiple shots)

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
	
	SGN_TS_PARS * fwdtspars = (SGN_TS_PARS *)(fwdfdm->fdpars);
	SGN_TS_PARS * bwdtspars = (SGN_TS_PARS *)(bwdfdm->fdpars);
	// copy dt and tspars from ref model to pert model
	fwdfdm->parcopy(bwdtspars,fwdtspars);
	(bwdmodel->tsind).dt = (fwdmodel->tsind).dt;
	
	// initializing  mig_model 
	int dim;         /* problem dimension */
	RDOM * ddom = &(bwdmodel->ld_a);
	rd_ndim(ddom,D_MP0,&dim);
	
	if (printact) fprintf(stream,"ASGAdjSampler::init - xrec=%d zeroing output\n",trace.t.tg.xrec) ;
	//        ra_a_zero(&(ddom->_s[D_MP0]));
	for (int i=0;i<dim;i++) 
	  ra_a_zero(&(ddom->_s[D_MV[i]]));

	IPNT ran;       
	size_t ntot;    
	
	rd_size(ddom,D_MP0,ran);
	if (printact) fprintf(stream,"ASGAdjSampler::init - D_MP0 dim=%d n1=%d n2=%d n3=%d\n",dim,ran[0],ran[1],ran[2]);
	ntot=1;
	for (int i=0;i<dim;i++) ntot*=ran[i];
	
	for (int j=0;j<(int)ntot;j++) ddom->_s[D_MP0]._s0[j]= REAL_ZERO;
	
	/*

	  for (int i=0;i<dim;i++) {
	  rd_size(ddom,D_MV[i],ran);
	  int m=1;
	  if ( err ) {
	  RVLException e;
	  e<<"Error: ASGAdjSampler::flush()\n";
	  e<<"Error: ra_setnull cannot apply on tmp array\n";
	  throw e;
	  }
	  for (int j=0;j<dim;j++) m*=ran[j];
	  for (int j=0;j<m;j++) ddom->_s[D_MV[i]]._s0[j]=REAL_ZERO;
	  }
	*/
      }// end of 'if(!err)...' branch 
    } /* end of static_init block */
    
    // fprintf(stream,"ASGAdjSampler::init, before sampler_init, err =%d\n",err);
    //    cerr<<"ASGAdjSampler::init, before sampler_init, err = "<<err<<endl;
    if (!err)
      err=sampler_init(&trace,
		       bwdmodel,
		       &(state.getPAR()),
		       state.getStream());

    // this should be the only terminal output in the entire run
    if (!err) cerr<<
                "adj ->"<<
		" rkw="<<retrieveGlobalRank()<<
		" rkl="<<retrieveRank()<<
		" src="<<trace.t.tg.irec<<endl;
    if(!err) is_samplerinit = true;
    
    //cerr<<"ASGAdjSampler::init, after sampler_init, err = "<<err<<endl;

    //fprintf(stream,"ASGAdjSampler::init, before aps.init, err =%d\n",err);
    fflush(stream);
    //cerr<<"ASGAdjSampler::init, before aps.init, err = "<<err<<endl;
    if (!err)
      aps.init();
   
    //cerr<<"ASGAdjSampler::init, after aps.init, err = "<<err<<endl;
    
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
      // fprintf(stderr,"in ASGAdjSampler::init: nt=%d istop=%d dt=%e\n",trace.t.tg.nt,trace.t.istop,trace.t.tg.dt);
      //  tt.setTargetTime(trace.t.istop,trace.t.tg.dt);
      //  bwdtt.setTargetTime(trace.t.istart,trace.t.tg.dt);
      TSIndex _tstart(t);
      TSIndex _tend(te);
      tt.setTargetTime(_tend);
      bwdtt.setTargetTime(_tstart);
    }
  
    //    fprintf(stream,"\n---exit ASGAdjSampler::init()\n ");
    //    fflush(stream);
   
    if (printact) fprintf(stream,"ASGAdjSampler::init, exit with err = %d\n",err);
     
    if (err) return false;
    return true;
  }

  void ASGAdjSampler::flush() const {

    // this is a no-op if model is physical (not extended)
    // and simulation is not finished
    // if physical model, invoke barrier to ensure that all contributions 
    // are available to reduction
    FILE *stream = state.getStream(); 
    //    fprintf(stream,"# bins = %d irec = %d nrecs = %d first = %d last = %d\n",get_panelnum_grid(state.getIWAVE().model.g),this->getRecIndx(),this->getNumRecs(),this->getFirstRec(),this->getLastRec());
    fflush(stream);
    //    cerr<<"# bins = "<<get_panelnum_grid(state.getIWAVE().model.g)<< " rec idx == "<<this->getRecIndx() <<" nrecs = "<<this->getNumRecs()<<"\n";
    if (get_panelnum_grid(state.getIWAVE().model.g) == 1) {
      //      if (this->getRecIndx() < this->getNumRecs()-1) return;
      if (this->getRecIndx() < this->getLastRec()) return;
#ifdef IWAVE_USE_MPI
      MPI_Barrier(retrieveGlobalComm());
#endif
    }
    PARARRAY & pars = state.getPAR();
    char * migbulkname = NULL;   // pert-bulkmodulus from migration
    char * migbuoyname = NULL;  // pert-density from migration
    char * ctmp = NULL;

    //    fprintf(state.getStream(),"---------ASGAdjSampler::flush, write out migrated models for data record = %d\n",this->getRecIndx());
       
    //ps_flcstring(pars,"mkappa",&migbulkname);
    //ps_flcstring(pars,"mbuoyancy",&migbuoyname);

    if(!ps_flcstring(pars,"mvelocity",&migbulkname)){
      ctmp = migbulkname;
      if(!ps_flcstring(pars,"mkappa",&migbulkname)) {
	userfree_(ctmp);
	ctmp = migbulkname;
	if(!ps_flcstring(pars,"mbulkmod",&migbulkname)){
	  userfree_(ctmp);
	  fprintf(state.getStream(),"\n NOTE: mkappa, mbulkmod and mvelocity are all supplied; mbulkmod is used \n");
	}
	else {
	  fprintf(state.getStream(),"\n NOTE: both mkappa and mvelocity are supplied; mkappa is used \n");
	}
      }
      else if(!ps_flcstring(pars,"mbulkmod",&migbulkname)){
	userfree_(ctmp);
        fprintf(state.getStream(),"\n NOTE: both mbulkmod and mvelocity are supplied; mbulkmod is used \n");
      }
      else {
	userfree_(migbulkname);
	//	cerr<<"\n ERROR: only find mvelocity; mbulkmod or mkappa need to be supplied \n";
	RVLException e;
	e<<"Error: ASGAdjSampler::flush()\n";
	e<<"ERROR: only find mvelocity; mbulkmod or mkappa need to be supplied \n";
	throw e;
      }
    }
    else if(!ps_flcstring(pars,"mkappa",&migbulkname)){
      ctmp = migbulkname;
      if(!ps_flcstring(pars,"mbulkmod",&migbulkname)){ 
	userfree_(ctmp);
	fprintf(state.getStream(),"\n NOTE: both mbulkmod and mkappa are supplied; mbulkmod is used \n");
      }
    }
    else 
      ps_flcstring(pars,"mbulkmod",&migbulkname);
 
    if(!ps_flcstring(pars,"mdensity",&migbuoyname)){
      ctmp = migbuoyname;
      if(!ps_flcstring(pars,"mbuoyancy",&migbuoyname)){
	userfree_(ctmp);
	fprintf(state.getStream(),"\n NOTE: both mdensity and mbuoyancy are supplied; mbuoyancy is used\n");
      }
      else {
	userfree_(migbuoyname);
	//	cerr<<"\n ERROR: only find mdensity; mbuoyancy need to be supplied \n";
	RVLException e;
	e<<"Error: ASGAdjSampler::flush()\n";
	e<<"ERROR: only find mdensity; mbuoyancy need to be supplied \n";
	throw e;
      }
    }
    else ps_flcstring(pars,"mbuoyancy",&migbuoyname);
     
    if (!migbulkname) {
      RVLException e;
      e<<"Error: ASGAdjSampler::flush()\n";
      e<<"failed to extract value for key=mbulkmod from param table\n";
      throw e;
    }
    if (!migbuoyname) {
      RVLException e;
      e<<"Error: ASGAdjSampler::flush()\n";
      e<<"failed to extract value for key=mbuoyancy from param table\n";
      throw e;
    }

    int err = 0;
    int dim;
    IPNT rags, rage;
    IPNT raVgs, raVge;
    IPNT ran;
        
    rd_ndim(&(state.getIWAVE().model.ld_a),D_MP0,&dim);
    rd_size(&(state.getIWAVE().model.ld_a),D_MP0,ran);
    rd_gse(&(state.getIWAVE().model.ld_a),D_MP0,rags,rage);

    // for non-extended models, reduce group adjoint accumulations to global
#ifdef IWAVE_USE_MPI
    /*
    cerr<<"----------ASGAdjSampler::flush bulkmod="<<migbulkname
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
      err=MPI_Reduce(state.getIWAVE().model.ld_a._s[D_MP0]._s0,buf,ntot,IWAVE_MPI_REAL,MPI_SUM,0,retrieveRemComm());
      if (err) {
	RVLException e;
	e<<"Error: AsgAdjSampler::flush() from MPI_Reduce, D_MP0\n";
	throw e;
      }
      // on group 0, copy back to alloca domain
      if (retrieveGroupID()==0) memcpy(state.getIWAVE().model.ld_a._s[D_MP0]._s0,buf,ntot*sizeof(ireal));
      delete [] buf;
    }
   
    int g1 = get_panelnum_grid(state.getIWAVE().model.g);
    int g2 = retrieveGroupID();
    // write occurs only if (i) extended modeling - all groups, or (ii) group 0, otherwise
    if ((g1 > 1) || (g2==0)) {

#endif

      //      cerr<<"---------ASGAdjSampler::flush - write bulkmod gid="<<retrieveGroupID()<<endl; 
      //      fprintf(state.getStream(),"---------ASGAdjSampler::flush migrated bulkmod = %s record=%d\n",migbulkname,this->getRecIndx());
      err = rsfwrite(state.getIWAVE().model.ld_a._s[D_MP0]._s0,
		     rags,
		     ran,
		     migbulkname,
		     0,
		     state.getStream(),
		     this->getRecIndx()); 
    
      //      cerr<<"-----------ASGAdjSampler::flush, wrote migrated bulk-modulus to "<<migbulkname<<", err = "<<err<<endl;
      if (err) {
	RVLException e;
	e<<"Error: ASGAdjSampler::flush()\n";
	e<<"rsfwrite return err for writing mig-dkappa\n";
	throw e;
      }

#ifdef IWAVE_USE_MPI
    }
#endif

    // for non-extended models, reduce group adjoint accumulations to global
#ifdef IWAVE_USE_MPI
    /*
      cerr<<"----------ASGAdjSampler::flush buoyancy="<<migbuoyname
      <<" rkw="<<retrieveGlobalRank()
      <<" rkl="<<retrieveRank()
      <<" gid="<<retrieveGroupID()
      <<endl;
    */
    if (!(get_panelnum_grid(state.getIWAVE().model.g) > 1)) {
      // foreach shifted buoyancy
      for (int i=0; i<dim;i++) {
	
	rd_size(&(state.getIWAVE().model.ld_a),D_MV[i],ran);
	// allocate buffer for reduction target
	int ntot=1;
	for (int j=0;j<dim;j++) ntot*=ran[j];
	ireal * buf = new ireal[ntot];
	
	err=MPI_Reduce(state.getIWAVE().model.ld_a._s[D_MV[i]]._s0,buf,ntot,IWAVE_MPI_REAL,MPI_SUM,0,retrieveRemComm());
	if (err) {
	  RVLException e;
	  e<<"Error: AsgAdjSampler::flush() from MPI_Reduce, D_MP["<<i<<"]\n";
	  throw e;
	}
	//	cerr<<"on group 0, copy back to alloca domain"<<endl;
	if (retrieveGroupID()==0) memcpy(state.getIWAVE().model.ld_a._s[D_MV[i]]._s0,buf,ntot*sizeof(ireal));
	delete [] buf;
	//	cerr<<"------------ axis "<<i<<"done"<<endl;
      }
    }

    // write occurs only if (i) extended modeling - all groups, or (ii) group 0, otherwise
    if ((get_panelnum_grid(state.getIWAVE().model.g) > 1) || retrieveGroupID()==0) {
#endif

      //    cerr<<"---------ASGAdjSampler::flush - write buoyancy gid="<<retrieveGroupID()<<endl; 
        
      IPNT ip;
      ireal q, qq;

      // build integer grid output array for buoyancy        
      RARR Btmp;
      IPNT ragsb, rageb, ranb;
      rd_gse(&(state.getIWAVE().model.ld_a), D_MP0, ragsb, rageb);
      /* pad bulkmod array to allow for adj-averaging shifted buoyancies */
      for (int i=0; i<dim; i++) { ragsb[i]--; rageb[i]++; }

      err = ra_create(&Btmp, dim, ragsb, rageb);
      if ( err ) {
	fprintf(state.getStream(),
		"Error: ASG_Sampler::flush - cannot allocate tmp array.\n");
	RVLException e;
	e<<"Error: ASG_Sampler::flush from ra_create\n";
	e<<"cannot allocate tmp array, err="<<err<<"\n";
	throw e;
      }
      ra_size(&Btmp, ranb);
      ra_zero(&Btmp);

      for (int i = 0; i < dim; i++) {
	rd_gse(&(state.getIWAVE().model.ld_a),D_MV[i],raVgs,raVge);

	for (int j = dim; j < RARR_MAX_NDIM; ++j ){
	  ragsb[j] = rageb[j] = 0;
	  rags[j] = rage[j] = 0;
	  raVgs[j] = raVge[j] = 0;
	} 
	/*
	  cerr<<"dim = "<<i<<endl;
	  cerr<<"raVgs =["<<raVgs[0]<<", "<<raVgs[1]<<", "<<raVgs[2]<<"]"<<endl;
	  cerr<<"raVge =["<<raVge[0]<<", "<<raVge[1]<<", "<<raVge[2]<<"]"<<endl;
	*/
#if RARR_MAX_NDIM > 2
	for ( ip[2] = raVgs[2] ; ip[2] <= raVge[2] ; ++ip[2] ) {
#endif
#if RARR_MAX_NDIM > 1
	  for ( ip[1] = raVgs[1] ; ip[1] <= raVge[1] ; ++ip[1] ) {
#endif
	    // adjoint-average D_MV[i] in direction i onto Btmp
	    for ( ip[0] = raVgs[0] ; ip[0] <= raVge[0] ; ++ip[0] ) {
	      //	    fprintf(stderr,"e: i=%d ip=[%d,%d] ->gget D_MV\n",i,ip[0],ip[1]);
	      //	    cerr<<"gget D_MV["<<i<<"]="<<D_MV[i]<<" ip=["<<ip[0]<<", "<<ip[1]<<"]\n";
	      q = ra_gget(&(state.getIWAVE().model.ld_a._s[D_MV[i]]), ip);
	      //	    cerr<<"gget Btmp ip=["<<ip[0]<<", "<<ip[1]<<"]\n";
	      qq = ra_gget(&Btmp, ip) + q * 0.5;
	      //	    cerr<<"gset Btmp ip=["<<ip[0]<<", "<<ip[1]<<"]\n";
	      ra_gset(&Btmp, ip, qq);
	      ip[i]++;
	      //	    cerr<<"gget Btmp ip=["<<ip[0]<<", "<<ip[1]<<"]\n";
	      qq = ra_gget(&Btmp, ip) + q * 0.5;
	      //	    cerr<<"gset Btmp ip=["<<ip[0]<<", "<<ip[1]<<"]\n";
	      ra_gset(&Btmp, ip, qq);
	      ip[i]--;
	    }
#if RARR_MAX_NDIM > 1
	  }
#endif
#if RARR_MAX_NDIM > 2
	}
#endif
      }

      //    cerr<<"\n---------- after averaging ----------------\n";

      if (state.getIWAVE().printact > 5) {
      
      	fprintf(state.getStream(),"\n----------------before rsfwrite-------------------\n");
	for (ip[1]=0;ip[1]<ranb[1];ip[1]++) {
	  for (ip[0]=0;ip[0]<ranb[0];ip[0]++) {
	    fprintf(state.getStream(),"%12.4e ",(Btmp._s0)
		    [ip[0]+ip[1]*ranb[0]]);
	  }
	  fprintf(state.getStream(),"\n");
	}
      }
      //      fprintf(state.getStream(),"---------ASGAdjSampler::flush migrated buoyancy = %s\n",migbuoyname);
      err = rsfwrite(Btmp._s0,
		     ragsb,
		     ranb,
		     migbuoyname,
		     0,
		     state.getStream(),
		     this->getRecIndx());
    
      //    cerr<<"ASGAdjSampler::flush, wrote migrated buoyancy to "<<migbuoyname<<", err = "<<err<<endl;
      if (err) {
	RVLException e;
	e<<"Error: ASGAdjSampler::flush()\n";
	e<<"rsfwrite return err "<<err<<" for writing mig-dbuoyancy\n";
	throw e;
      }

      //      fprintf(state.getStream(),"\n----------------returned from rsfwrite-------------------\n");
      // deallocate buoyancy integer grid buffer
      ra_destroy(&Btmp);

      // matches isextmodel or group ID = 0
#ifdef IWAVE_USE_MPI
    }
#endif

    /* reset is_samplerinit to false for next record */
    is_samplerinit = false;

    userfree_(migbulkname);
    userfree_(migbuoyname);

    //    fprintf(state.getStream(),"-----------ASGAdjSampler::flush exit\n");
    //    fflush(state.getStream());
  } 


  ostream & ASGAdjSampler::write(ostream & str) const {
    sampler_fprint(&trace,state.getStream());
    aps.fprint();
    return str;
  }

}

