#include "istate.hh"

#define IWAVE_VERBOSE

namespace TSOpt {

  void IWaveEnvironment(int argc, char ** argv, int ts,
			PARARRAY ** par, 
			FILE ** stream) {
    int err=0;
    //    cerr<<"IWaveEnvironment -> initparallel_global\n";
    initparallel_global(ts);
    //    cerr<<"IWaveEnvironment -> initinoutstream\n";
    if ((err=initoutstream(stream,retrieveGlobalRank(),retrieveGlobalSize()))) {
      RVLException e;
      e<<"Error: IWaveEnvironment from initoutstream, err="<<err<<"\n";
      e<<"failed to initialize output stream\n";
      throw e;
    }
    //    cerr<<"IWaveEnvironment -> readinput\n";
#ifdef IWAVE_VERBOSE
    for (int i=1;i<argc;i++) 
      cerr<<"IWaveEnvironment: argv["<<i<<"]="<<argv[i]<<endl;
#endif
    if ((err=readinput(par,*stream,argc,argv))) {
      RVLException e;
      e<<"Error: IWaveEnvironment from readinput, err="<<err<<"\n";
      e<<"failed to parse param file \n";
      throw e;
    }
#ifdef IWAVE_VERBOSE
    cerr<<"IWaveEnvironment: after construct par\n";
    ps_printall(*(*par),stderr);
#endif
    //    cerr<<"IWaveEnvironment -> initparallel_local\n";
    if ((err=initparallel_local(*(*par),*stream))) {
      RVLException e;
      e<<"Error: IWaveEnvironment from initparallel_local, err="<<err<<"\n";
      e<<"failed to create cartesian grid or local comm\n";
      throw e;
    }
    //    cerr<<"IWaveEnvironment -> exit\n";
  }

  // helper function - not accessible outside this file!
  void synch(IWAVE * pstate, 
	     bool fwd,
	     int it,
	     int iv,
	     IWaveInfo const & ic,
	     FILE * stream) {
    
    // no-op unless MPI is defined
    
#ifdef IWAVE_USE_MPI
    
    int err = 0;
    // step info for the state to be updated

    int ia;       /* array index */

    int i;
    MPI_Status status;  
    int tmpdest, tmpsource;
    MPI_Datatype tmpdest_dt, tmpsource_dt;
    int tmpdest_val, tmpsource_val;
    void *tmpdest_buf, *tmpsource_buf;
    double time;
    time = 0.0; /* To avoid "unitialized variable" warning */

    fprintf(stream,"IWAVE::synch: printact = %d\n",pstate->printact);

    /* we use the time step internal index for the perturbed field because
       this is the field being updated - on call should be same as index for
       unperturbed field, but don't use this */
    if ( (pstate->printact > 5) ) {
      fprintf(stream,"\n------ giwave_synch: before exchange, step %d substep %d\n",it,iv);
      for (ia=0;ia<RDOM_MAX_NARR;ia++) {
	if (fd_update(ia,iv,ic)) {
	  fprintf(stream,"------ iarr = %d\n",ia);
	  rd_print(&((pstate->model).ld_a), ia, stream);
	}
      }
      fflush(stream); 
    }

    for (ia=0;ia<RDOM_MAX_NARR;ia++) {
      if (fd_update(ia,iv,ic)) {
	if ( (pstate->printact > 1) ) {
	  fprintf(stream,"\n------ giwave_synch fwd=%d array=%d -------------\n",fwd,ia);
	  fflush(stream); 
	}

	IPNT gs;
	IPNT ge;
	RARR rsave;
	ra_setnull(&rsave);

	for ( i = 0; i < (pstate->model).nnei; ++i ) {
	  // create send data - corresponds to ld_s	
	  if ( (pstate->pinfo).seinfo[ia][i].type != MPI_DATATYPE_NULL ) {  
	    tmpdest = (pstate->pinfo).sranks[i];
	    tmpdest_buf = (pstate->pinfo).seinfo[ia][i].buf;
	    tmpdest_dt = (pstate->pinfo).seinfo[ia][i].type;
	    // adjoint case - save a copy of send RARR, which becomes recv buffer
	    // for adjoint

	    if (!fwd) {
	      rd_gse(&(((pstate->model).ld_s)[i]),ia,gs,ge);
	      if ((err=ra_create(&rsave,gs,ge))) {
		fprintf(stream,"\nError: giwave_synch from ra_create err=%d\n",err);
		RVLException e;
		e<<"\nError: IwaveSynch from ra_create err="<<err<<"\n";
		throw e;
	      }
	      ra_zero(&rsave);
	      if ((err=ra_copy(&rsave,&(((((pstate->model).ld_s)[i])._s)[ia])))) {
		fprintf(stream,"\nError: giwave_synch from ra_copy err=%d\n",err);
		RVLException e;
		e<<"\nError: IWaveSynch from ra_copy err="<<err<<"\n";
		throw e;
	      }
	    }

	  }	  
	  else {
	    tmpdest = MPI_PROC_NULL;
	    tmpdest_buf = &tmpdest_val;
	    tmpdest_dt = MPI_INT;
	  } 
	  // create receive data - corresponds to ld_r
	  if ( (pstate->pinfo).reinfo[ia][i].type != MPI_DATATYPE_NULL ) {
	    tmpsource = (pstate->pinfo).rranks[i];
	    tmpsource_buf = (pstate->pinfo).reinfo[ia][i].buf;
	    tmpsource_dt = (pstate->pinfo).reinfo[ia][i].type;
	  } 
	  else {
	    tmpsource = MPI_PROC_NULL;
	    tmpsource_buf = &tmpsource_val;
	    tmpsource_dt = MPI_INT;
	  }

	  if ((pstate->printact > 1)) {
	    fprintf(stream, "    i = %d, sending to wrk = %d, receiving from wrk = %d, [NULL = %d]\n", 
		    i, tmpdest, tmpsource, MPI_PROC_NULL); 
	    fflush(stream);
	  }
	  cerr<<"s/r\n";	
	  /* 
	     "dest" is send buffer
	     "source" is receive buffer

	     if fwd: 
	     SEND data in tmpdest_buf to rk=tmpdest
	     RECEIVE data in tmpsource_buf from rk=tmpsource

	     else:
	     SEND data in tmpsource_buf to rk=tmpsource
	     RECEIVE data in tmpdest_buf from rk=tmpdest

	  */
	  cerr<<"s/r fwd="<<fwd<<"\n";
	  if ((tmpdest != MPI_PROC_NULL) && (tmpsource != MPI_PROC_NULL)) {
	    if (fwd) {
	      err = MPI_Sendrecv(tmpdest_buf, 1, tmpdest_dt, tmpdest, iv,
				 tmpsource_buf, 1, tmpsource_dt, tmpsource, iv,
				 (pstate->pinfo).ccomm, &status);
	      if ( err != MPI_SUCCESS ) {
		fprintf(stream, 
			"ERROR. Internal: MPI_Sendrecv error #%d, nei=%d, iv=%d, arr=%d. ABORT.\n", 
			err, i, iv, ia);
		RVLException e;
		e<<"    ERROR. Internal: MPI_Sendrecv error #"<<err<<", nei="<<i<<", iv="<<iv<<", arr="<<ia<<" ABORT\n";
		throw e;	      
	      }
	    }
	    else {
	      // first receive into send buffer, which has been copied into a tmp buffer
	      // reciprocally, send receive buffer
	      err = MPI_Sendrecv(  tmpsource_buf, 1, tmpsource_dt, tmpsource, iv,
				   tmpdest_buf, 1, tmpdest_dt, tmpdest, iv,
				   (pstate->pinfo).ccomm, &status);
	      
	      if ( err != MPI_SUCCESS ) {
		fprintf(stream, 
			"ERROR. Internal: MPI_Sendrecv error #%d, nei=%d, iv=%d, arr=%d. ABORT.\n", 
			err, i, iv, ia);
		RVLException e;
		e<<"    ERROR. Internal: MPI_Sendrecv error #"<<err<<", nei="<<i<<", iv="<<iv<<", arr="<<ia<<" ABORT\n";
		throw e;	      
	      }

	      // then add tmp buffer back
	      if ( ( (pstate->pinfo).seinfo[ia][i].type != MPI_DATATYPE_NULL ) ) {
		if (!(rsave._s0)) {
		  fprintf(stream,"\nError: giwave_synch before axpy: rsave not initialized\n");
		  RVLException e;
		  e<<"\nError: IWaveSynch before axpy: rsave not initialized\n";
		  throw e;
		}
		fflush(stream);
		if ((err=ra_axpy(&(((((pstate->model).ld_s)[i])._s)[ia]),&rsave,REAL_ONE))) {
		  fprintf(stream,"\nError: giwave_synch from ra_axpy err=%d\n",err);
		  ra_dump(&(((((pstate->model).ld_r)[i])._s)[ia]),stream);
		  ra_dump(&rsave,stream);
		  ra_destroy(&rsave);
		  RVLException e;
		  e<<"Error: IWaveSynch from ra_axpy err="<<err<<"\n";
		  throw e;
		}
		ra_destroy(&rsave);
	      }

	    // now comp subdom data is correct - synchronize
	    /////// 13.04.14: apparently unnecessary, but zeroing the receive buffer IS!!
	    /*
	    err = MPI_Sendrecv(tmpdest_buf, 1, tmpdest_dt, tmpdest, iv,
			       tmpsource_buf, 1, tmpsource_dt, tmpsource, iv,
			       (pstate->pinfo).ccomm, &status);
	    if ( err != MPI_SUCCESS ) {
	      fprintf(stream, 
		      "ERROR. Internal: MPI_Sendrecv error #%d, nei=%d, iv=%d, arr=%d. ABORT.\n", 
		      err, i, iv, ia);
	      RVLException e;
	      e<<"    ERROR. Internal: MPI_Sendrecv error #"<<err<<", nei="<<i<<", iv="<<iv<<", arr="<<ia<<" ABORT\n";
	      throw e;
	    }
	    */
	    // now zero the receive buffer
	      if ( ( (pstate->pinfo).reinfo[ia][i].type != MPI_DATATYPE_NULL ) ) { 
		ra_zero(&(((((pstate->model).ld_r)[i])._s)[ia])); 
	      }
	    }
	  }
	}
      }
    }

    if ( (pstate->printact > 5) ) {
      fprintf(stream,"\n------ giwave_synch: after exchange, step %d substep %d\n",it,iv);
      for (ia=0;ia<RDOM_MAX_NARR;ia++) {
	if (fd_update(ia,iv,ic)) {
	  fprintf(stream,"------ iarr = %d\n",ia);
	  rd_print(&((pstate->model).ld_a), ia, stream);
	}
      }
      fflush(stream); 
    }
#endif
  }

  IWaveSampler::IWaveSampler(IWAVE * state, string key, PARARRAY & pars, FILE * stream)
    : axes(0), prev_panelindex(-1), tg(NULL), dump_term(0) {
    
    samplekey=key;
    pname="";
    if (parse(pars,key,pname)) {
      string gname = "";

      // extract proto name - use preferentially to actual file name,
      // particularly to handle temp filenames, which do not have suffixes.
      // note that only geometrical info will be extracted here, which can
      // as well be found from the proto files
      const char * rname = NULL;

      // because current implementation of out-of-core classes does all 
      // file handling on rk 0, must extract prototype file name there.
      // this block of code should be abtracted eventually to hide 
      // implementation details against future changes eg. to distrtibuted
      // or network database. As is, will work properly so long as prototype
      // file is used in iwave_fopen on rk=0.

#ifdef IWAVE_USE_MPI
      int cprotolen=0;
      if (retrieveGlobalRank()==0) {
#endif
	if ((rname = iwave_getproto(pname.c_str()))) gname = rname;
	else gname = pname;
#ifdef IWAVE_USE_MPI
	// C string length 
	cprotolen = gname.size()+1;
      }

      MPI_Bcast(&cprotolen,1,MPI_INT,0,retrieveGlobalComm());
      char * cproto = new char[cprotolen];
      if (retrieveGlobalRank()==0) strcpy(cproto,gname.c_str());
      MPI_Bcast(cproto,cprotolen,MPI_CHAR,0,retrieveGlobalComm());
      gname=cproto;
      delete [] cproto;
#endif

      // check for file type
      size_t pos = gname.find_last_of(".");
      if (pos==std::string::npos || pos >= gname.size()-1) {
	RVLException e;
	e<<"Error: IWaveSampler constructor - filename "<<gname<<" has no suffix\n";
	e<<"gname = "<<gname<<"\n";
	e<<"pname = "<<pname<<"\n";
	e<<"rname = "<<rname<<"\n";
	throw e;
      }
      size_t net = gname.size()-pos-1;
      suffix = gname.substr(pos+1,net);
#ifdef IWAVE_VERBOSE      
      cerr<<"IWaveSampler: gname="<<gname<<" pname="<<pname<<" rank="<<retrieveGlobalRank()<<endl;
#endif

      if (suffix=="rsf") {
	grid g;
	if (read_grid(&g,gname.c_str(),stream)) {
	  RVLException e;
	  e<<"Error: IWaveSampler constructor - rsf case, failed to read grid\n";
	  e<<"  key = "<<samplekey<<" proto filename = "<<gname<<"\n";
	  throw e;
	}
	// copy in rigid order - 0 - dim-1 are spatial axes in order (z, x, y); 
	// dim is time, hence skipped; all others are extended - identify by axis id
	// also check that the index mapping is bijective
	for (int i=0; i<g.gdim; i++) {
	  for (int j=0;j<i;j++) {
	    if (g.axes[i].id == g.axes[j].id) {
	      RVLException e;
	      e<<"Error: IWaveSampler constructor\n";
	      e<<"  repetition of axis ids\n";
	      throw e;
	    }
	  }
	}
	for (int k=0; k<g.gdim; k++) {
	  // don't include internal extended axes! they do not
	  // contribute to simulation grid - they are ** internal **
	  if (g.axes[k].id < EXTINT) {
	    //	    cerr<<"  copy axis "<<i<<"\n";
	    axis * a = new axis;
	    copy_axis(a,&(g.axes[k]));
	    axes.push_back(a);
	  }
	  else {
	    //	    cerr<<"  do not copy axis "<<i<<"\n";
	  }
	}

	// determine whether the grid has spatial axes - these would have id's 
	// between 0 and spatial dim
	has_Spatial_Axes = true;
	for (int j=0; j<(state->model).g.dim; j++) 
	  has_Spatial_Axes = has_Spatial_Axes && has_Axis(j);

	//      cerr<<"rsf case: axis vector = \n";
	//      for (int i=0;i<axes.size();i++) {
	//	fprint_axis(stderr,*(axes[i]));
	//      }
	//      cerr<<"exit IWaveSampler constructor, rsf case\n";
      }
      else if (suffix=="su") {

	// by default, not used to sample spatial fields - 
	// only trace time series
	has_Spatial_Axes = false;

	//      cerr<<"in IWaveSampler constructor, su case\n";

	/* construct trace geometry object */
	tg = new tracegeom;
	int err=construct_tracegeom(tg,
				    pname.c_str(),
				    (state->model).tsind.dt,
				    SRC_TOL,
				    stream);
  
	if (err) {
	  fprintf(stream,"ERROR: IWaveSampler\n");
	  fprintf(stream,"return from construct_tracegeom with err=%d\n",err);
	  RVLException e;
	  e<<"ERROR: IWaveSampler\n";
	  e<<"  return from construct_tracegeom with err="<<err<<"\n";
	  throw e;
	}

	dump_term = valparse<int>(pars,"dump_term",0);
    
	/* set up axes */
	axis * at = new axis;
	at->n = tg->nt;
	at->o = tg->t0;
	at->d = (state->model).tsind.dt;
	at->id = (state->model).g.dim;
	axes.push_back(at);

	/*
	  if (retrieveGlobalRank()==1) {
	  cerr<<"IWaveSampler - su case, time axis: n="<<at->n<<" d="<<at->d<<" o="<<at->o<<" id="<<at->id<<endl;
	  cerr<<"IWaveSampler - su case, model grid =\n";
	  fprint_grid(stderr,(state->model).g);
	  cerr<<"IWaveSampler - su case, new axis =\n";
	  fprint_axis(stderr,*at);
	  cerr<<"IWaveSampler - su case, recorded axis "<<axes.size()-1<<"=\n";
	  fprint_axis(stderr,*(axes[axes.size()-1]));
	  cerr<<"IWaveSampler - su case, nrec ="<<tg->nrec<<"\n";
	  }
	*/

	tracestart=(int)(floor(((at->o)/(at->d))+0.1));
	tracestop=tracestart + at->n - 1;

	if (tg->nrec > 1) {
	  axis * ae = new axis;
	  ae->n = tg->nrec;
	  ae->o = 0.0;
	  ae->d = 1.0;
	  ae->id = (state->model).g.dim + 1;
	  axes.push_back(ae);
	}
	
	/* read sampling order */
	sampord=0;
	parse(pars,"sampord",sampord);
	
      }
      else {
	RVLException e;
	e<<"Error: IWaveSampler constructor - suffix "<<suffix<<" represents unknown data structure\n";
	e<<"  rk="<<retrieveGlobalRank()<<"\n";
	iwave_fprintall(stderr);
	throw e;
      }
    }

  }

  IWaveSampler::~IWaveSampler() {
    if (tg) {
      destroy_tracegeom(tg);
      delete tg;
    }
    for (int i=0;i<(int)axes.size(); i++) {
      delete axes[i];
    }
  }

  axis const & IWaveSampler::getAxis(int i) const { 
    try {
      return *(axes.at(i)); 
    }
    catch (out_of_range) {
      RVLException e;
      e<<"Error: IWaveSampler::getAxis\n";
      e<<"  axis index "<<i<<" out of range = [0, "<<axes.size()<<"\n";
      throw e;
    }
  }

  bool IWaveSampler::has_Axis(int id) const {
    bool has = false;
    for (int i=0; i<(int)axes.size(); i++) 
      has = has || (axes[i]->id == id);
    return has;
  }

  axis const * IWaveSampler::getTimeAxis(int dim) const {
    for (int i=0; i<(int)axes.size();i++) {
      if (axes[i]->id == dim) return axes[i];
    }
    return NULL;
  }

  ireal IWaveSampler::getCellVol() const {
    ireal cv = REAL_ONE;
    for (int i=0; i<(int)axes.size();i++) 
      cv *= axes[i]->d;
    return cv;
  }

  ireal IWaveSampler::getRecipCellVol() const {
    ireal cv = REAL_ONE;
    ireal rcv = REAL_ONE;
    for (int i=0; i<(int)axes.size();i++) 
      cv *= axes[i]->d;
    if (ProtectedDivision<ireal>(REAL_ONE,cv,rcv)) {
      RVLException e;
      e<<"Error: IWaveSampler::getRecipCellVol\n";
      e<<"  samplekey = "<<samplekey<<"\n";
      e<<"  zerodivide by cell volume - check data geometry\n";
      throw e;
    } 
    return rcv;
  }
    
  void IWaveSampler::sample(grid g, IPNT step, bool fwd, bool input, 
			    IWAVE * state, int ridx, int iwdx, FILE * stream,
			    bool dryrun, ostream & drystr) {
    try {

      if (this->axes.size() > 0) {

	bool load = false;       // load data from disk
	bool save = false;       // save data to disk
	bool init = false;       // initialize internal buffers and transfer info, if needed
	bool sample = false;     // sample to internal or external buffer

	// compute istart
	IPNT istart;
	IPNT istop;
	get_gs(istart,g);
	get_ge(istop,g);

	// detect time axis
	axis const * a = this->getTimeAxis(g.dim);

	// compute panelindex = step through external extended indices
	// load - if at BEGIN of all axes that you DON'T HAVE, and at begin of time axis
	// save - if at END of all axes that you DON'T HAVE, and at end of time axis
	// for determining "begin" and "end" of time axis, sense is reversed when
	// fwd=false

	int panelnum=1;
	int panelindex=0;

	// nontrivial extended axes
	if (g.gdim > g.dim+1) {
	  // multipanel movies not permitted
	  /*
	  if (has_Spatial_Axes && a) {
	    RVLException e;
	    e<<"Error: IWaveSampler::sample, samplekey="<<samplekey<<"\n";
	    e<<"  time axis i/o of full spatial grid\n";
	    e<<"  not allowed in conjunction with multiple\n";
	    e<<"  simulations (no multishot movies)\n";
	    throw e;
	  }
	  */
	  for (int i=g.dim+1; i<g.gdim; i++) {
	    panelindex += (step[i]-istart[i])*panelnum;
	    panelnum *= istop[i]-istart[i]+1;
	    //	    cerr<<"axis "<<i<<" has="<<has_Axis(g.axes[i].id) << " input="<<input<<"\n";
	  }
	}
	else {
	  panelnum=1;
	  panelindex=0;
	}

	// compute group panelindex limits
	int first;
	int last;
	calc_group(&first, &last, panelnum);
	if ((panelindex < first) ||
	    (panelindex > last)) {
	  RVLException e;
	  e<<"Error: IWaveSampler::sample\n";
	  e<<"  panelindex "<<panelindex<<" out of range ["<<first<<", "<<last<<"]\n";
	  throw e;
	}
#ifdef IWAVE_VERBOSE
	cerr<<"IWaveSampler::sample: samplekey="<<samplekey<<"\n";
	cerr<<"  first="<<first<<" last="<<last<<"\n";
	cerr<<"  panelindex="<<panelindex<<"\n";
#endif

	// now compute sampler-specific panel index

	int this_panelnum=1;
	int this_panelindex=prev_panelindex;	

	// case 1: extended axes, excluding time axis
	if (dryrun)
	  drystr<<"\nSampler decision loop: samplekey="<<samplekey<<"\n";
	std::vector<int> myaxes(0);
	for (int i=g.dim+1; i<g.gdim; i++) {
	  if (has_Axis(g.axes[i].id)) {
	    myaxes.push_back(i);
	  }
	}
	if (myaxes.size()>0) this_panelindex=0;
	for (size_t i=0; i< myaxes.size(); i++) {
	  this_panelindex += (step[myaxes[i]]-istart[myaxes[i]])*this_panelnum;
	  this_panelnum *= istop[myaxes[i]]-istart[myaxes[i]]+1;
	}
	if (dryrun) 
	  drystr<<myaxes.size()<<" extended axes, sampling enabled\n";
	// no extended axes
	if (myaxes.size()==0) {
	  if ((input && (panelindex==first)) ||
	      (!input && (panelindex==last))) {
	    this_panelindex++;
	    this_panelindex=iwave_min(this_panelindex, this_panelnum-1);
	    if (dryrun) 
	      drystr<<"no extended axes - sampling enabled, update this_panelindex"<<endl;
	  }
	  else {
	    this_panelindex=prev_panelindex;
	    if (dryrun)
	      drystr<<"no extended axes - sampling disabled, do not update this_panelindex"<<endl;
	  }	  
	}

	if (a && has_Spatial_Axes) {
	  // now account for time panels, if time axis is part of this sampler's
	  // grid and it samples a full spatial field - then sampling and load/save
	  // are assumed to be the same
	  this_panelnum = a->n;
	  // compute timeslices for sampling
	  float t = g.axes[g.dim].d*(step[g.dim]-istart[g.dim]) + g.axes[g.dim].o;
	  if (dryrun) drystr<<"  sampler on rarr "<<ridx<<" iwave "<<iwdx<<" time = "<<t<<"\n";
	  // warn if sampling will be forced to duplicate
	  if ((istop[g.dim]-istart[g.dim]+1 < (int)(a->n)) && (step[g.dim]==istart[g.dim])) {
	    cerr<<"NOTE: IWaveSampler::sample\n";
	    cerr<<"  number of time steps = "<<istop[g.dim]-istart[g.dim]+1;
	    cerr<<" less than requested number of time samples = "<<a->n<<"\n";
	    cerr<<"  some time samples will be duplicated\n";
	  }
	  // time before begin of sample range
	  if (t < a->o - 0.01*a->d || t > a->o + (a->n-1+0.01)*(a->d)) this_panelindex = prev_panelindex;
	  else {
	    // compute relative time on axis with positive margin of one-half internal time step
	    float curr_time_rel   = ((t - a->o + 0.5*g.axes[g.dim].d)/(a->d));
	    if (dryrun) drystr<<"curr_time_rel="<<(int)(floor(curr_time_rel))<<"\n";
	    // get panel index
	    this_panelindex = iwave_max(0,iwave_min((int)(floor(curr_time_rel)),this_panelnum-1));
	  }
	}	    

	// time condition for load/save
	// note that input flags input for forward mode - has
	// opposite meaning for reverse mode
	bool buffered = (a && !has_Spatial_Axes);
	bool atbeg = (fwd && (step[g.dim]==istart[g.dim])) ||
	  (!fwd && (step[g.dim]==istop[g.dim]));
	bool atend = (!fwd && (step[g.dim]==istart[g.dim])) ||
	  (fwd && (step[g.dim]==istop[g.dim]));

	if (dryrun && (this_panelindex==prev_panelindex)) {
	  drystr<<"this_panelindex unchanged - no sampling:\n";
	  drystr<<"this_panelindex = "<<this_panelindex<<" this_panelnum="<<this_panelnum<<endl; 
	}
	if (this_panelindex != prev_panelindex) {

	  if (dryrun) {
	    drystr<<"this_panelindex has changed - sampling is possible\n";
	    drystr<<"this_panelnum="<<this_panelnum<<" this_panelindex="<<this_panelindex<<" prev_panelindex="<<prev_panelindex<<"\n";
	    drystr<<"input="<<input<<" buffered="<<buffered<<" atbeg="<<atbeg<<" atend="<<atend<<"\n";
	  }

	  if (buffered) {
	    // for buffered i/o, load only at start, save only at 
	    // end of time axis
	    if (atbeg && input) {
	      load=true;
	      prev_panelindex = this_panelindex;
	    }	   
	    if (atend && !input) {
	      save=true;
	      prev_panelindex = this_panelindex;
	    }
	  }
	  // for unbuffered i/o, load/save whenever panelindex changes, if
	  // you have a time axis, otherwise only at ends of time axis
	  // also zero target array at ends of time axis, not otherwise
	  else if (a) {
	    load = input;
	    save = !input;
	    prev_panelindex = this_panelindex;	    
	  }
	  else {
	    load = input && atbeg;
	    save = !input && atend;
	    if (load || save)  
	      prev_panelindex = this_panelindex;
	    if (load) {
	      if (int err = ra_a_zero(&((state->model).ld_a._s[ridx]))) {
		RVLException e;
		e<<"Error: IWaveSampler::sample from ra_a_zero\n";
		e<<"err="<<err<<"\n";
		throw e;
	      }
	    }
	  }
	}

	// for buffered i/o, sample at every time step, 
	// initialize buffer at begin of time axis
	if (buffered) {
	  sample=true;
	  init = atbeg;
	}
	// otherwise, sampling is synonymous with load/save
	else {
	  sample = (input && load) || (!input && save);
	}

	// internal grid info
	int dim;
	IPNT n0;
	IPNT gs0;
	IPNT ge0;
	IPNT n;
	IPNT gs;
	IPNT ge;
	RDOM * dom = &((state->model).ld_a) ;
	if (rd_ndim(dom,ridx,&dim) ||
	    rd_size(dom,ridx,n0) ||
	    rd_gse(dom,ridx,gs0,ge0)) {
	  RVLException e;
	  e<<"Error: IWaveSample::sample\n";
	  e<<"  failure on size query, allocated array\n";
	  e<<"  iwdx="<<iwdx<<" ridx="<<ridx<<"\n";
	  throw e;
	}
	if (rd_size(&((state->model).ld_c),ridx,n) ||
	    rd_gse(&((state->model).ld_c),ridx,gs,ge)) {
	  RVLException e;
	  e<<"Error: IWaveSample::sample\n";
	  e<<"  failure on size query, computational array\n";
	  e<<"  iwdx="<<iwdx<<" ridx="<<ridx<<"\n";
	  throw e;
	}

	// data buffer for current RARR
	float * data = dom->_s[ridx]._s0;
	  
	// scale for adjoint linearized i/o
	//   - input - by cell volume 
	//   - output - by reciprocal cell volume
	// detect linearized i/o by iwdx - if > 0 this i/o op
	// reads or writes a perturbation
	ireal scale = REAL_ONE;
	if (!fwd && input && (iwdx>0)) scale = this->getCellVol();
	if (!fwd && !input && (iwdx>0)) scale = this->getRecipCellVol();
	
	if (sample) {

	  if (dryrun) {
	    
	    drystr<<"\nIWaveSample: samplekey="<<samplekey<<"\n";
	    drystr<<"  panelnum        = "<<panelnum<<"\n";
	    drystr<<"  panelindex      = "<<panelindex<<"\n";
	    drystr<<"  this_panelnum   = "<<this_panelnum<<"\n";
	    drystr<<"  this_panelindex = "<<this_panelindex<<"\n";
	    drystr<<"  input           = "<<input<<"\n";
	    drystr<<"  fwd             = "<<fwd<<"\n";
	    drystr<<"  load            = "<<load<<"\n";
	    drystr<<"  save            = "<<save<<"\n";
	    drystr<<"  init            = "<<init<<"\n";
	    if (fwd) {
	      drystr<<"  fwd time step    = "<<step[g.dim]<<"\n";
	    }
	    else {
	      drystr<<"  bwd time step    = "<<step[g.dim]<<"\n";
	    }
	    if (step[g.dim]==istart[g.dim]) drystr<<"  at start of time loop\n";
	    if (step[g.dim]==istop[g.dim])  drystr<<"  at end of time loop\n";
	    if (load) drystr<<"  read data\n";
	    if (save) drystr<<"  write data\n";
	    for (int i=g.dim+1;i<g.gdim;i++) {
	      if (this->has_Axis(g.axes[i].id)) {
		drystr<<"  extd axis    = "<<i<<" id = "<<g.axes[i].id<<" step = "<<step[i]<<"\n";
	      }
	    }
	    
	    if (this_panelnum > 1) {
	      drystr<<"  this_panelnum     = "<<this_panelnum<<"\n";
	      drystr<<"  this_panelindex   = "<<this_panelindex<<"\n";
	    }
	    drystr<<"  iwave index  = "<<iwdx<<"\n";
	    drystr<<"  rarray index = "<<ridx<<"\n";
	    if (a) drystr <<"  has time axis\n";
	    if (has_Spatial_Axes) drystr<<"  samples spatial fields\n";
	  }
	  
	  else {
	    if (suffix=="rsf") {
	      if (load) {
#ifdef IWAVE_VERBOSE
		cerr<<"rsf load on samplekey="<<this->samplekey<<" iwdx="<<iwdx<<" ridx="<<ridx<<"\n";
		cerr<<"sample 48 before read = "<<data[48]<<endl;
#endif
		// 03.02.14: read into buffer; axpy via loop - hoist all scaling
		// out of rsfread, which is used in overwrite mode so no need to
		// independently initialize buf

		int err=rsfread(data,gs0,n0,pname.c_str(),scale,stream,this_panelindex);
		if (err) {
		  RVLException e;
		  e<<"Error: IWaveSample::sample, rsf case\n";
		  e<<"  failed to read panel "<<panelindex<<" from "<<pname<<"\n";
		  throw e;
		}	      
		// 03.02.14: temporarily disable scaling for unit conversion - scale
		// only for adjoint input, update array by saxpy
#ifdef IWAVE_VERBOSE
		cerr<<"sample 48 after read = "<<data[48]<<endl;
#endif	      
	      }  
	      if (save) {
#ifdef IWAVE_VERBOSE
		cerr<<"rsf save on samplekey="<<this->samplekey<<" iwdx="<<iwdx<<" ridx="<<ridx<<endl;
#endif
		int err=0;

		// stack case 
		if (this_panelnum == 1) {

#ifdef IWAVE_USE_MPI
		  // if no time axis, then local stack
		  MPI_Barrier(MPI_COMM_WORLD);
		  size_t ntot=1;
		  for (int i=0;i<dim;i++) ntot *= n0[i];
		  float * buf = new float[ntot];
		  // reduce along remote comm axis - same domain in every group
		  int err=MPI_Reduce(data,buf,ntot,IWAVE_MPI_REAL,MPI_SUM,0,retrieveRemComm());
		  if (err) {
		    RVLException e;
		    e<<"Error: IWaveSampler::sample from MPI_Reduce\n";
		    e<<"  iwdx = "<<iwdx<<" ridx = "<<ridx<<"\n";
		    throw e;
		  }
		  // on group 0, copy back to alloca domain
		  if (retrieveGroupID()==0) {
		    memcpy(data,buf,ntot*sizeof(ireal));
		  }
		  delete [] buf;
		  if (retrieveGroupID()==0) {
#endif
		  err=rsfwrite(data,gs0,n0,pname.c_str(),scale,stream,0);
		  if (err) {
		    RVLException e;
		    e<<"Error: IWaveSample::sample, rsf single panel / time sample case\n";
		    e<<"  failed to write panel global index = "<<panelindex<<"\n";
		    e<<"  local index = "<<this_panelindex<<" from "<<pname<<"\n";
		    e<<"  return value from rsfwrite = "<<err<<"\n";
		    throw e;
		  }
#ifdef IWAVE_USE_MPI
		  }
#endif
		}
		else {
		  err=rsfwrite(data,gs0,n0,pname.c_str(),scale,stream,this_panelindex);
		  if (err) {
		    RVLException e;
		    e<<"Error: IWaveSample::sample, rsf multiple panel case\n";
		    e<<"  failed to write panel global index = "<<panelindex<<"\n";
		    e<<"  local index = "<<this_panelindex<<" from "<<pname<<"\n";
		    e<<"  return value from rsfwrite = "<<err<<"\n";
		    throw e;
		  }
		  // if more than one panel index for this sampler, then
		  // zero out accumulation in connected rarray, in preparation
		  // for computing the next panel
		  // probably suffices for most instances
		  // BUT don't do this if sampling on the time axis!!!
		  if (!a) ra_a_zero(&(dom->_s[ridx]));
		}
	      }
	    }

	    else if (suffix=="su") {
	      
	      /* extract grid params - space and time */
	      IPNT nl;                       /* axis lengths, local grid */
	      RPNT ol;                       /* axis origins, local grid */
	      RPNT dl;                       /* axis steps, local grid */
	      IPNT ng;                       /* axis lengths, global grid */
	      RPNT og;                       /* axis origins, global grid */
	      RPNT dg;                       /* axis steps, global grid */
	      IPNT axord;                    /* axis order array */
	      
	      get_n(nl,(state->model).gl);
	      get_o(ol,(state->model).gl);
	      get_d(dl,(state->model).gl);
	      get_n(ng,(state->model).g);
	      get_o(og,(state->model).g);
	      get_d(dg,(state->model).g);
	      get_ord(axord,(state->model).g);

	      // int flag for initializing trace buffer
	      // note that "load" is a subcase of "init" for buffered i/o
	      int initbuf = 0;
	      if (!fwd && load) initbuf = -1; // -1 -> adjoint cubic interpolation
	      if ( fwd && load) initbuf =  1; //  1 -> cubic interpolation

#ifdef IWAVE_VERBOSE
	      if (load) cerr<<"su load on samplekey="<<samplekey<<endl;
	      if (save) cerr<<"su save on samplekey="<<samplekey<<endl;
#endif

	      if (init) {
		//		cerr<<"samplekey="<<samplekey<<" rk="<<retrieveGlobalRank()<<" init_tracegeom\n";
		/* initialize tracegeom - note data is UNSCALED */
		int err=init_tracegeom(tg,
				       this_panelindex,
				       og,nl,dl,ol,axord,
				       sampord,
				       (state->model).g.dim,
				       initbuf,
				       stream);
		if (err) {
		  RVLException e;
		  e<<"Error: IWaveSampler::sample from init_tracegeom\n";
		  e<<"  err = "<<err<<"\n";
		  throw e;
		}

		if (dump_term) fprint_tracegeom(tg,stream);

	      }

	      // sample if in axis range
	      if ((step[g.dim] >= tracestart) &&
		  (step[g.dim] <= tracestop)) {

		// read/write trace switch		
		int traceinput = 0;
		if (input) {
		  traceinput = 1;
		  // this assumes that for reference state
		  // that (i) all trace input is RHS (questionable),
		  // and (ii) all trace input should be treated as
		  // defining a numerical delta (safe)
		  if (iwdx==0) {
		    ireal cv=1.0f;
		    for (int kk=0;kk<(state->model).g.dim;kk++) 
		      cv*=(state->model).g.axes[kk].d;
		    scale = (state->model).tsind.rhs/cv;
		    //		    cerr<<"rhs="<<(state->model).tsind.rhs<<" cv="<<cv<<" iwdx="<<iwdx<<" scale="<<scale<<endl;
		  }
		}

		sampletraces(tg,
			     sampord,
			     traceinput,
			     step[g.dim] - tracestart,
			     n0, gs0,
			     n, gs,
			     data,
			     scale);
	      }

	      if (save) {
		//		cerr<<"samplekey="<<samplekey<<" rk="<<retrieveGlobalRank()<<" writetraces\n";
		int err=writetraces(tg,dg,og,stream);
		if (err) {
		  RVLException e;
		  e<<"Error: IWaveSampler::sample from writetraces\n";
		  e<<"  err = "<<err<<"\n";
		}
	      }
	    }
	  }
	}

      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from IWaveSampler::sample\n";
      throw e;
    }
  }
    
  IWaveTree::IWaveTree(std::vector<IWAVE *> sv, IWaveInfo const & _ic)
    : own(false), sa(iwave_max(sv.size()/2,0)),
      rd(iwave_max(sv.size()/2,0)), ic(_ic) {
    // cerr<<"iwavetree private constructor\n";
    if (sv.size()>1) {
      for (size_t i=0;i<sv.size()/2;i++) {
	sa[i]=sv[i];
	rd[i]=&((sv[i]->model).ld_c);
      }
    }
    else {
      RVLException e;
      e<<"Error: IWaveTree private constructor\n";
      e<<"  input vector of size "<<sv.size()<<"\n";
      throw e;
    }
    ref=NULL;
    // cerr<<"exit iwavetree private constructor\n";
  }

  IWaveTree::IWaveTree(PARARRAY & _pars, FILE * _stream, IWaveInfo const & _ic,
		       int order)
    : own(true), sa(pow2(order)), rd(pow2(order)), ic(_ic) {
    try {
      // cerr<<"iwavetree constructor\n";
      for (size_t i=0; i<sa.size(); i++) {
	sa[i]=new IWAVE;
	int err=0;
	// cerr<<"i="<<i<<" call iwave_construct\n";
	if ((err=iwave_construct(sa[i],&(_pars),_stream,ic))) {
	  RVLException e;
	  e<<"Error: IWaveTree main constructor, IWAVE["<<i<<"]\n";
	  e<<"  returning error "<<err<<" from iwave_construct\n";
	  throw e;
	}
	// clean it up
	// cerr<<"-> rd_a_zero\n";
	if ((err=rd_a_zero(&((sa[i]->model).ld_a)))) {
	  RVLException e;
	  e<<"Error: IWaveTree main constructor, IWAVE["<<i<<"]\n";
	  e<<"  returning error "<<err<<" from rd_a_zero\n";
	  throw e;
	}
	rd[i]=&((sa[i]->model).ld_c);
      }
      // create reference state out of first half of
      // array, without allocating any new memory,
      // using private constructor
      // cerr<<"priv constr\n";
      if (order>0) ref = new IWaveTree(sa,ic);
      else ref=NULL;
      // extract time step, set up step index 
      // cerr<<"exit IWaveTree constructor\n";
    }
    catch (RVLException & e) {
      e<<"\ncalled from IWaveTree constructor\n";
      throw e;
    }
  }

  IWaveTree::~IWaveTree() {
    if (ref) delete ref;
    if (own) 
      for (size_t i=0; i<sa.size(); i++) {
	iwave_destroy(sa[i],ic.get_mdest());
	delete sa[i];
      }    
  }
   
  IWaveSim::IWaveSim(int _order, 
		     bool _fwd, 
		     PARARRAY & pars, 
		     FILE * _stream, 
		     IWaveInfo const & _ic, 
		     int _printact,
		     int _snaps,
		     bool _dryrun,
		     ostream & _drystr,
		     ostream & _announce)
    : ic(_ic), fwd(_fwd), stream(_stream),  
      printact(_printact), order(_order), snaps(_snaps),
      cps(NULL), narr(0), 
      dryrun(_dryrun), drystr(_drystr), 
      announce(_announce) {
    try {
      // cerr<<"iwavesim constr\n";
      // sanity
      if (!fwd && ((snaps<=0) || (order<=0))) {
	RVLException e;
	e<<"Error: IWaveSim constructor\n";
	e<<"  adjoint mode (fwd=false):"<<"\n";
	if (snaps<=0) {
	  e<<"  must provide positive number of \n";
	  e<<"  storage units for checkpoints (argument \"snaps\")\n";
	  e<<"  however snaps = "<<snaps<<"\n";
	}
	if (order<=0) {
	  e<<"  does not make sense for reference simulation (order=0)\n";
	  e<<"  order supplied in argument list = "<<order<<"\n";
	}
	throw e;
      }

      
      // cerr<<"step 1: create list of i/o tasks\n";
      IOTask(t,order,fwd,ic);
#ifdef IWAVE_VERBOSE
      cerr<<"IWaveSim constructor: fwd="<<fwd<<" order="<<order<<endl;
      IOTaskWriter(t,cerr);
#endif
      // cerr<<"step 2: build state\n";
      w = new IWaveTree(pars, stream, ic, order);

      // cerr<<"step 2a: build grid.\n";
      // start with spatial grid, which has been initialized
      // in particular g.dim = problem spatial dimn
      // copy only the spatial axes, leaving the rest to
      // the sampler constructors
      init_default_grid(&g);
      int dim = (((w->getStateArray())[0])->model).g.dim;      
      int check = 0;
      for (int i=0;i<RARR_MAX_NDIM;i++) {
	int id = (((w->getStateArray())[0])->model).g.axes[i].id; 
	if ( (id> -1) && (id < dim) ) {
	  copy_axis(&(g.axes[id]),&((((w->getStateArray())[0])->model).g.axes[i]));
	  check++;
	}
      }
      if (check != dim) {
	RVLException e;
	e<<"Error: IWaveSim constructor\n";
	e<<"  grid in model indices funged up\n";
	throw e;
      }
      // trash all axes above spatial dimn
      for (int i=dim;i<RARR_MAX_NDIM;i++) {
	g.axes[i].n=0;
	g.axes[i].d=0.0;
	g.axes[i].o=0.0;
	g.axes[i].id=-1;
      }
      g.dim =  dim;
      g.gdim=g.dim;
      // cerr<<"---------------\n";
      // cerr<<"IWaveSim: sim grid g after initial construction from model.g:\n";
      //      fprint_grid(stderr,g);
      // cerr<<"---------------\n";

      // axis order is rigidly:
      // z, x, y (for 3d) or z, x (for 2d);
      // axis index g.dim is time
      // axis indices above g.dim are extended
      // note that this ordering has nothing to do with actual 
      // storage order in rsf data structure - correspondence is
      // by axis id
      
      // cerr<<"step 2a: build checkpoint structure if required\n";
      // work out number of dynamic arrays
      ndyn = 0;
      for (int i=0;i<RDOM_MAX_NARR;i++) {
	// cerr<<"i="<<i<<" ndyn="<<ndyn<<endl;
	if (fd_isarr(i,(w->getStateArray())[0]->model,ic) && fd_isdyn(i,ic)) ndyn++; 
      }
      if (ndyn==0) {
	RVLException e;
	e<<"Error: IWaveSim constructor\n";
	e<<"  count of dynamic fields = 0\n";
	throw e;
      }

      // cerr<<"step 2b: in adjoint case, allocate checkpoint buffers\n";
      // IWAVEs 0,...,order-1 contain the reference data for the 
      // adjoint computation, so need order * snaps * ndyn RARRAYs.
      if (!fwd) {
	typedef RARR* RP;
	typedef RARR** RPP;
	narr = snaps*pow2(order-1)*ndyn;
	RARR * cpstmp = new RARR[narr];
	RARR ** cpstmp2 = new RP[snaps*pow2(order-1)];
	for (size_t i=0;i<snaps*pow2(order-1);i++) cpstmp2[i]=&(cpstmp[i*ndyn]);
	cps = new RPP[snaps];
	for (int i=0;i<snaps;i++) cps[i]=&(cpstmp2[i*pow2(order-1)]);
	int l=0;
	for (int k=0;k<(w->getRefStateArray())[0]->model.ld_a.narr;k++) {
	  if (fd_isarr(k,(w->getRefStateArray())[0]->model,ic) && fd_isdyn(k,ic)) {
	    // pull out gs, ge for kth rarr in w->getRefStateArray[0] 
	    IPNT gs;
	    IPNT ge;
	    ra_a_gse(&((w->getRefStateArray())[0]->model.ld_a._s[k]),gs,ge);
	    //	    int ndim = (w->getRefStateArray())[0]->model.ld_a._s[k].ndim;
	    for (int i=0;i<snaps;i++) {
	      for (size_t j=0;j<pow2(order-1);j++) {
		// ra_create cps[i][j][l] -- check l
		if (l<0 || l>ndyn-1) {
		  RVLException e;
		  e<<"Error: IWaveSim: construct checkpoint array\n";
		  e<<"  bad news - index into rarray "<<l<<" out of range [0,"<<ndyn<<"\n";
		  throw e;
		}
		if (int err=ra_create(&(cps[i][j][l]),gs,ge)) {
		  RVLException e;
		  e<<"Error: IWaveSim from ra_create - checkpoint array\n";
		  e<<"  snap # = "<<i<<"\n";
		  e<<"  branch # = "<<j<<"\n";
		  e<<"  dyn arr # = "<<l<<"\n";
		  e<<"  rarr # = "<<k<<"\n";
		  e<<"  gs: ";
		  for (int i=0; i<RARR_MAX_NDIM; i++) e<<gs[i]<<" ";
		  e<<"\n";
	       	  e<<"  ge: ";
		  for (int i=0; i<RARR_MAX_NDIM; i++) e<<ge[i]<<" ";
		  e<<"\n";
		  e<<"  err = "<<err<<"\n";
		  throw e;
		}
	      }
	    }
	    l++;
	  }
	}
      }

      // cerr<<"step 3: construct list of samplers, axes\n";
      // note sampler exposed to first IWAVE - grid, time step info
      s.clear();
      for (size_t i=0; i<t.size(); i++) {
	IWaveSampler * tmp = NULL;
#ifdef IWAVE_VERBOSE
	cerr<<"construct sampler "<<i<<" on keyword "<<t[i]->keyword<<"\n";
#endif
	tmp = new IWaveSampler(w->getStateArray()[0], t[i]->keyword, pars, stream);
	// mod of 07.12.13: if current job does not provide source/sink for 
	// this i/o task, then sampler returns no axes. in that case, set
	// this pointer in s to NULL
	if (tmp->getNumAxes() == 0) {
#ifdef IWAVE_VERBOSE
	  cerr<<"sampler "<<i<<" on keyword "<<t[i]->keyword<<" is null\n";
#endif
	  delete tmp;
	  tmp = NULL;
	}
	else {
	  /////// loop through sampler axes, add each to grid, then on
	  /////// success add sampler address to sim list
	  //	  cerr<<"sampler "<<i<<" on keyword "<<t[i]->keyword<<" is live\n";
	  // note these are internal simulation axes, not archival, 
	  // which will require some re-arranging in traceio
	  // grid_union(grid *, std::vector<axis *>)
	  //	  cerr<<"in sampler loop: sampler index = "<<i<<" keyword = "<<t[i]->keyword<<"\n";
	  //	  for (int k=0;k<tmp->getNumAxes(); k++) 
	  //	    fprint_axis(stderr,tmp->getAxis(k));
	  //	  fprintf(stderr,"  add to grid:\n");
	  //	  fprint_grid(stderr,g);

	  // IMPORTANT MOD 11.12.14: since physical domain is spec'd by FIELD ID 0,
	  // only add spatial axes for FIELD ID 0. These axes were already fixed in 
	  // step 1, by initial construction of grid. So skip all subsequent axes
	  // with id < dim

	  for (int j=0; j<tmp->getNumAxes();j++) {
	    // IMPORTANT CHANGE 07.12.13: the time axis is special, because
	    // the internal step is already fixed by the IWAVE constructor.
	    // therefore, change this axis to have the simulation time step
	    // before registering it via grid_union.
	    // Change 08.12.13: copy axis to workspace - only const access
	    // now allowed
	    axis * a = new axis;
	    init_default_axis(a);
	    copy_axis(a, &(tmp->getAxis(j)));
	    float tntmp = 0.0f;
	    if (a->id == g.dim) {     // signals time axis!!
	      if (ProtectedDivision<float>((a->n * a->d),
					   (w->getStateArray()[0]->model.tsind.dt),
					   tntmp)) {
		RVLException e;
		e<<"Error: IWaveSim constructor\n";
		e<<"  zerodivide by dt stored in IWaveTree model struct\n";
		e<<"  dt = "<<w->getStateArray()[0]->model.tsind.dt<<"\n";
		throw e;
	      }
	      // adjust number, spacing
	      a->n = iwave_max(1, (int)(tntmp+0.1f));
	      a->d = w->getStateArray()[0]->model.tsind.dt;
	      // adjust origin so that in all cases 0 is a (virtual) grid point
	      // this assures that all time axes are compatible. Adjust o and n
	      // so that new grid interval contains old.
	      float tmpo = a->o;
	      a->o =
		//		((int)(((a->o)/(a->d))-0.1))
		((int)((a->o)/(a->d)))
		*(a->d);
	      if (tmpo<a->o) {
		a->o -= a->d;
		a->n ++;
	      }
	      // notice that now simulation time axis has been DETERMINED, so 
	      // must be used in trace sampler!!!!
	    }

	    if ((a->id > g.dim-1) && (!grid_union(&g,a))) {
	      //	      cerr<<"Error: IWaveSim constructor from grid_union\n";
	      RVLException e;
	      e<<"Error: IWaveSim constructor from grid_union\n";
	      fprintf(stream,"Error: IWaveSim constructor from grid_union\n");
	      fprintf(stream,"  failed to add these axes:\n");
	      for (int k=0;k<tmp->getNumAxes(); k++) 
		fprint_axis(stream,tmp->getAxis(k));
	      fprintf(stream,"  to grid:\n");
	      fprint_grid(stream,g);
	      throw e;
	    }
	    // having added to grid, or not, trash it
	    delete a;
	  
	  }
	  //      	  fprintf(stderr,"  result grid:\n");
	  //       	  fprint_grid(stderr,g);
	}
	s.push_back(tmp);
      }
      // 08.01.14: panelindex construction makes this unnecessary
      //      fprint_grid(stderr,g);
      // step 4: modify grid if necessary
      // no-op if no mpi
      //#ifdef IWAVE_USE_MPI
      //      mpi_update_grid(&g);
      //#endif
      // cerr<<"exit iwavesim constructor\n";
    }
    catch (RVLException & e) {
      e<<"\ncalled from IWaveSim constructor\n";
      throw e;
    }
  }
  
  IWaveSim::~IWaveSim() {
    //    cerr<<"destructor\n";
    for (size_t i=0; i<t.size(); i++) {
      //cerr<<"destroy sampler "<<i<<endl;
      if (s.at(i)) {
	//	cerr<<"destroy sampler "<<i<<endl;
	delete s.at(i);
      }
    }
    //    cerr<<"destroy state\n";
    if (w) delete w;
    for (size_t i=0; i<t.size(); i++) {
      //      cerr<<"destroy iotask "<<i<<endl;
      if (t.at(i)) delete t.at(i); 
    }
    if (cps) {
      //      cerr<<"destroy checkpoints\n";
      for (int i=0;i<snaps;i++) {
	for (size_t j=0;j<pow2(order-1);j++) {
	  for (int k=0;k<ndyn;k++) {
	    ra_destroy(&(cps[i][j][k]));
	  }
	}
      }
      delete [] &(cps[0][0][0]);
      delete [] &(cps[0][0]);
      delete [] &(cps[0]);
    }
  }

  void IWaveSim::run() {
    try {
      // cerr<<"iwavesim::run\n";
      // initialize step
      IPNT step;
      IPNT start;
      IPNT stop;
      get_gs(start, g);
      get_ge(stop, g);
      // initially, assign step = start
      IASN(step,start);
      // next, work out first and last panels for extended
      // axes, this process
      int first;
      int last;
      int panelnum=1;
      int panelindex=0;
      for (int i=g.dim+1; i<g.gdim; i++) {
	panelnum *= stop[i]-start[i]+1;
      }
      calc_group(&first, &last, panelnum);
      // increment step until at first panel
      while (panelindex < first) {
	// update step array
	for (int i=g.dim+1;i<g.gdim;i++) {
	  if (step[i]<stop[i]) step[i]++;
	  else {
	    if (i<g.gdim-1) step[i]=start[i];
	  }
	}
	panelindex++;
      }

      // pull out fdpars for use in time step - same in every
      // step, and for every RDOM, so do it once here and get 
      // from root RDOM
      void * fdm = NULL;
      if (!dryrun) {
	fdm = (w->getStateArray()[0]->model.specs);
      }

      // time axis - index=dim: last data point if fwd = false (reverse)
      
      // step 5: initial sample - includes any post-contruction
      // initialization
      
      // step 6: sim loop - includes time, record,...
      // time step function knows what to do based 
      // merely on size. use next_step to increment
      // extended axes - time step loop is explicit

      // here more = more records
      for (panelindex=first; panelindex<=last; panelindex++) {
	/*
	for (int i=0;i<w->getStateArray().size();i++) {
	  if (int err=rd_a_zero(&((w->getStateArray()[i]->model).ld_a))) {
	    RVLException e;
	    e<<"Error: IWaveTree main constructor, IWAVE["<<i<<"]\n";
	    e<<"  returning error "<<err<<" from rd_a_zero\n";
	    throw e;
	  }
	}
	*/
	// compute panelindex = step through external extended indices
	if (dryrun) {
	  if (panelindex==first) {
	    drystr<<"\n*** IWaveSim: model="<<ic.iwave_model<<" fwd="<<fwd<<" deriv="<<order<<"\n";
	  }
	  drystr<<"    panel="<<panelindex<<" in range=["<<first<<", "<<last<<"], rank="<<retrieveGlobalRank()<<endl;
	}
	else {
	  if (panelindex==first && retrieveGlobalRank()==0) {
	    announce<<"\n*** IWaveSim: model="<<ic.iwave_model<<" fwd="<<fwd<<" deriv="<<order<<"\n";
	  }
	  announce<<"    panel="<<panelindex<<" in range=["<<first<<", "<<last<<"], rank="<<retrieveGlobalRank()<<endl;

#ifdef IWAVE_VERBOSE
	  announce<<"simulation grid:\n";
	  for (int i=0;i<g.gdim;i++) {
	    if (i<=g.dim) step[i]=start[i];
	    announce<<"axis "<<i<<" n="<<g.axes[i].n<<" o="<<g.axes[i].o<<" d="<<g.axes[i].d<<" id="<<g.axes[i].id<<" gs="<<start[i]<<" ge="<<stop[i]<<" step="<<step[i]<<"\n";
	  }
#endif
	}
	if (dryrun) {
	  drystr<<"\nIWaveSim::run - initialize dynamic arrays"<<endl;
	  drystr<<"simulation grid = \n";
	  for (int i=0;i<g.gdim;i++) {
	    if (i<=g.dim) step[i]=start[i];
	    drystr<<"axis "<<i<<" n="<<g.axes[i].n<<" o="<<g.axes[i].o<<" d="<<g.axes[i].d<<" id="<<g.axes[i].id<<" gs="<<start[i]<<" ge="<<stop[i]<<" step="<<step[i]<<"\n";
	  }
	  drystr<<"\n\n";
	}

	if (fwd) {
	  // cerr<<"\nIWaveSim::run - FORWARD WET RUN\n\n";
	  if (dryrun) {
	    drystr<<"\nIWaveSim::run - FORWARD DRY RUN\n\n";
	  }
#ifdef IWAVE_VERBOSE
	  cerr<<"IWaveSim::run - dynamic init\n";
#endif
	  for (size_t i=0;i<w->getStateArray().size();i++) {
	    iwave_dynamic_init(w->getStateArray()[i],start[g.dim],ic);
	  }

	  for (int it=start[g.dim]; it<stop[g.dim]; it++) {

#ifdef IWAVE_VERBOSE
	    cerr<<"IWaveSim::run - it="<<it<<"\n";
#endif
	    if (dryrun) drystr<<"\n";
	    step[g.dim]=it;

	    if (!dryrun && (printact > 0)) {
	      float dt = g.axes[g.dim].d;
	      float ot = g.axes[g.dim].o;
	      cerr<<"it="<<it<<" t="<<ot+it*dt<<endl;
	    }

	    
	    for (size_t i=0; i<t.size(); i++) {
	      if (s[i]) {
#ifdef IWAVE_VERBOSE
		cerr<<"IWaveSim::run - rk="<<retrieveGlobalRank()<<" sampler["<<i<<"]\n";
#endif
		s[i]->sample(g,step,fwd,
			     t[i]->input,w->getStateArray()[t[i]->iwaveindex],
			     t[i]->rarrindex,t[i]->iwaveindex,stream,
			     dryrun,drystr);
	      }
	    }
	    
	    // first time through, call check
#ifdef IWAVE_VERBOSE
	    cerr<<"IWaveSim::run - call check\n";
#endif
            if (!dryrun) {
	      if (it==start[g.dim]) ic.get_check()(w->getRDOMArray()[0],fdm,stream);
            }
	    if (dryrun) {
	      drystr<<"\nIWaveSim::run fwd step "<<it<<" -> "<<it+1<<"\n";
	    }
	    else { 
#ifdef IWAVE_VERBOSE
	      cerr<<"IWaveSim::run - rk="<<retrieveGlobalRank()<<": timestep\n";
#endif
	      for (int iv=0;iv<fd_numsubsteps(ic);iv++) {
		//		cerr<<"timestep\n";
		ic.get_timestep()(w->getRDOMArray(),fwd,iv,fdm);
		// in fwd loop, synch ALL dynamic arrays
		
		for (size_t k=0; k<w->getStateArray().size(); k++) {
		  //		  if (w->getStateArray()[k]->printact > 5) 
		  //		    fprintf(stream,"\n*** SYNCH: iwdx = %d\n\n",k);
#ifdef IWAVE_VERBOSE
		  cerr<<"IWaveSim::run -synch substep="<<iv<<" iwdx="<<k<<"\n";
#endif
		  synch(w->getStateArray()[k],fwd,it,iv,ic,stream);
		  //		  cerr<<"  k="<<k<<"\n";
		}
	      }
	    }
	  }
	  
	  if (dryrun) drystr<<"\n";
	  step[g.dim]=stop[g.dim];
	  for (size_t i=0; i<t.size(); i++) {
	    if (s[i]) {
#ifdef IWAVE_VERBOSE
	      cerr<<"IWaveSim::run - rk="<<retrieveGlobalRank()<<" sampler["<<i<<"]\n";
#endif
	      s[i]->sample(g,step,fwd,t[i]->input,
			   w->getStateArray()[t[i]->iwaveindex],
			   t[i]->rarrindex,t[i]->iwaveindex,stream,
			   dryrun,drystr);
	    }
	  }
	}

	else if (stop[g.dim]-start[g.dim]==1) {
	  
	  for (size_t i=0;i<w->getStateArray().size();i++) {
	    iwave_dynamic_init(w->getStateArray()[i],start[g.dim],ic);
	  }

	  int it = start[g.dim]; // ref state index
	  int at = stop[g.dim];  // pert state index
	  
	  if (dryrun) {
	    drystr<<"\nIWaveSim::run - ADJOINT DRY RUN\n";
	    drystr<<"  SINGLE STEP CASE\n\n";
	  }	 
#ifdef IWAVE_VERBOSE
	  cerr<<"\nIWaveSim::run - ADJOINT DRY RUN\n";
	  cerr<<"  SINGLE STEP CASE\n\n";
#endif
	  // load reference data
	  step[g.dim]=it;
	  bool reffwd = true;
	  for (size_t i=0; i<t.size(); i++) {
	    // sample data for reference
	    if (t[i]->iwaveindex < (int)pow2(order-1) && s[i]) {
	      s[i]->sample(g,step,reffwd,t[i]->input,
			   w->getStateArray()[t[i]->iwaveindex],
			   t[i]->rarrindex,t[i]->iwaveindex,stream,
			   dryrun,drystr);
	    }
	  }
	  // check
          if (!dryrun) {
	    ic.get_check()(w->getRDOMArray()[0],fdm,stream);
          }

	  // load adjoint data 
	  step[g.dim]=at;
	  for (size_t i=0; i<t.size(); i++) {
	    // pert sample
	    if (t[i]->iwaveindex >= (int)pow2(order-1) && s[i]) { 
	      s[i]->sample(g,step,fwd,t[i]->input,
			   w->getStateArray()[t[i]->iwaveindex],
			   t[i]->rarrindex,t[i]->iwaveindex,stream,
			   dryrun,drystr);
	    }
	  }
	  if (!dryrun) {
	    // backwards step - only need to synch top-order pert arrays
	    for (int iv=fd_numsubsteps(ic)-1; iv>-1; iv--) {
	      ic.get_timestep()(w->getRDOMArray(),fwd,iv,fdm);
	      for (size_t k=pow2(order-1); k<pow2(order); k++) {
		//		if (w->getStateArray()[k]->printact > 5) 
		//		  fprintf(stream,"\n*** SYNCH: iwdx = %d\n\n",k);
		synch(w->getStateArray()[k],fwd,at,iv,ic,stream);
	      }
	    }	
	  }
	  else {
	    drystr<<"IWaveSim::run adj - first adjoint step: at="<<at;
	  }
	  at--;
	  if (dryrun) {
	    drystr<<"->"<<at<<"; ref step = "<<it<<endl;
	  }
	  step[g.dim]=at;

	  for (size_t i=0; i<t.size(); i++) {
	    // pert sample
	    if (t[i]->iwaveindex >= (int)pow2(order-1) && s[i]) { 
	      s[i]->sample(g,step,fwd,t[i]->input,
			   w->getStateArray()[t[i]->iwaveindex],
			   t[i]->rarrindex,t[i]->iwaveindex,stream,
			   dryrun,drystr);
	    }
	  }

	}

	else {
	  // flag for step, sample on fwd steps
	  bool reffwd = true;

	  // create checkpoint schedule
	  ostringstream ostr;
	  Revolve r(stop[g.dim]-start[g.dim],snaps,ostr);

	  for (int i=0;i<snaps;i++) {
	    for (int j=0;j<(int)pow2(order-1);j++) {
	      for (int l=0;l<ndyn;l++) {
		if (int err=ra_a_zero(&(cps[i][j][l]))) {
		  RVLException e;
		  e<<"Error: IWaveSim from ra_create - zero rarray\n"; 
		  e<<"  snapindex="<<i<<" iwaveindex="<<j<<" rarrindex="<<l<<"\n";
		  e<<"  err = "<<err<<"\n";
		}
	      }
	    }
	  }

	  std::vector<int> cplist(snaps);
	  int it = start[g.dim]; // ref state index
	  int at = stop[g.dim];  // pert state index

#ifdef IWAVE_VERBOSE
	  cerr<<"\nIWaveSim::run adj - load initial data time step "<<it<<endl;
#endif
	  step[g.dim]=it;

	  for (size_t i=0;i<w->getStateArray().size();i++) {
	    iwave_dynamic_init(w->getStateArray()[i],start[g.dim],ic);
	  }
	  
	  for (size_t i=0; i<t.size(); i++) {
	    // sample data for reference
	    if (t[i]->iwaveindex < (int)pow2(order-1) && s[i]) {
	      s[i]->sample(g,step,reffwd,t[i]->input,
			   w->getStateArray()[t[i]->iwaveindex],
			   t[i]->rarrindex,t[i]->iwaveindex,stream,
			   dryrun,drystr);
	    }
	  }

	  // check if it=start
          if (!dryrun) {
	  if (it==start[g.dim]) ic.get_check()(w->getRDOMArray()[0],fdm,stream);
          }

	  ACTION::action whatodo;

	  if (dryrun) {
	    drystr<<"\nIWaveSim::run - ADJOINT DRY RUN\n";
	    drystr<<"  FORWARD TIME LOOP\n\n";
	  }
	  do {
	    whatodo = r.revolve();

	    if (whatodo == ACTION::takeshot) {
	      int cp = r.getcheck();
	      for (int j=0;j<(int)pow2(order-1);j++) {
		int l = 0;
		for (int k=0;k<RDOM_MAX_NARR;k++) {
		  if (fd_isarr(k,w->getStateArray()[0]->model,ic) && fd_isdyn(k,ic)) {		  
		    if (ra_a_copy(&(cps[cp][j][l]),&(((w->getRefRDOMArray())[j])->_s[k]))) {
		      RVLException e;
		      e<<"Error: IWaveSim::run\n";
		      e<<"attempt to store checkpoint "<<cp
		       <<" encountered error from ra_a_copy\n";
		      e<<"IWAVE index = "<<j<<" checkpoint RARR index = "
		       <<l<<" RDOM RARR index = "<<k<<"\n";
		      throw e;
		    }
		    l++;
		  }
		}
	      }
	      cplist.at(cp)=it;
	      if (dryrun) {
		drystr<<"\nIWaveSim::run adj - stored step "<<it<<" in checkpoint "<<cp<<endl;
	      }
#ifdef IWAVE_VERBOSE
	      cerr<<"\nIWaveSim::run adj - stored step "<<it<<" in checkpoint "<<cp<<endl;
#endif
	    }
	    if (whatodo == ACTION::advance) { 
	      for (int j = r.getoldcapo(); j < r.getcapo(); j++) {
		if (dryrun) drystr<<"\n";
		if (!dryrun) {
		  for (int iv=0;iv<fd_numsubsteps(ic);iv++) {
		    ic.get_timestep()(w->getRefRDOMArray(),reffwd,iv,fdm);
		    for (size_t k=0; k<pow2(order-1); k++) {
		      //		      if (w->getStateArray()[k]->printact > 5) 
		      //			fprintf(stream,"\n*** SYNCH: iwdx = %d\n\n",k);
		      synch(w->getRefStateArray()[k],reffwd,it,iv,ic,stream);
		    }
		  }
#ifdef IWAVE_VERBOSE	
		  cerr<<"IWaveSim::run adj - fwd step "<<it;
#endif
		}
		else {
		  drystr<<"IWaveSim::run adj - fwd step "<<it;
		}
		it++;
#ifdef IWAVE_VERBOSE
		cerr<<"->"<<it<<endl;
#endif
		if (dryrun) {
		  drystr<<"->"<<it<<endl;
		}
		step[g.dim]=it;
		for (size_t i=0; i<t.size(); i++) {
		  // sample data for reference
		  if (t[i]->iwaveindex < (int)pow2(order-1) && s[i]) {
		    s[i]->sample(g,step,reffwd,t[i]->input,
				 w->getStateArray()[t[i]->iwaveindex],
				 t[i]->rarrindex,t[i]->iwaveindex,stream,
				 dryrun,drystr);
		  }
		}
	      }
	    }
	    if (whatodo == ACTION::firsturn) {
#ifdef IWAVE_VERBOSE
	      cerr<<"\nIWaveSim::run - ADJOINT WET RUN\n";
	      cerr<<"  REVERSE TIME LOOP\n\n";
#endif
	      if (dryrun) {
		drystr<<"\nIWaveSim::run - ADJOINT DRY RUN\n";
		drystr<<"  REVERSE TIME LOOP\n\n";
	      }	      
	      step[g.dim]=at;
	      for (size_t i=0; i<t.size(); i++) {
		// need to clean output arrays at outset of reverse time loop
		if (t[i]->iwaveindex >= (int)pow2(order-1)) {
		  // iwave_dynamic_init(w->getStateArray()[t[i]->iwaveindex],at,ic);
		  //		  if (!(t[i]->input) )
		  //		    ra_a_zero(&(w->getStateArray()[t[i]->iwaveindex]->model.ld_a._s[t[i]->rarrindex]));
		}
		// pert sample 
		if (t[i]->iwaveindex >= (int)pow2(order-1) && s[i]) { 
		  s[i]->sample(g,step,fwd,t[i]->input,
			       w->getStateArray()[t[i]->iwaveindex],
			       t[i]->rarrindex,t[i]->iwaveindex,stream,
			       dryrun,drystr);
		}
	      }
	      if (!dryrun) {
		for (int iv=fd_numsubsteps(ic)-1; iv>-1; iv--) {
		  ic.get_timestep()(w->getRDOMArray(),fwd,iv,fdm);
		  for (size_t k=pow2(order-1); k<pow2(order); k++) {
		    //		    if (w->getStateArray()[k]->printact > 5) 
		    //		      fprintf(stream,"\n*** SYNCH: iwdx = %d\n\n",k);
		    synch(w->getStateArray()[k],fwd,at,iv,ic,stream);
		  }
		}	
#ifdef IWAVE_VERBOSE
		cerr<<"IWaveSim::run adj - first adjoint step: at="<<at;
#endif
	      }
	      else {
		drystr<<"IWaveSim::run adj - first adjoint step: at="<<at;
	      }
	      at--;
#ifdef IWAVE_VERBOSE
	      cerr<<"->"<<at<<"; ref step = "<<it<<endl;
#endif
	      if (dryrun) {
		drystr<<"->"<<at<<"; ref step = "<<it<<endl;
	      }
	    }
	    if (whatodo == ACTION::youturn) {
#ifdef IWAVE_VERBOSE
	      cerr<<"\n";
#endif
	      if (dryrun) drystr<<"\n";
	      step[g.dim]=at;
#ifdef IWAVE_VERBOSE
	      cerr<<"backwards step"<<at<<"before sample iwdx=1 ridx=1 ucb="<<
		(w->getRDOMArray()[order]->_s)[1]._s0[48]<<"\n";
#endif
	      for (size_t i=0; i<t.size(); i++) {
		// pert sample
		if (t[i]->iwaveindex >= (int)pow2(order-1) && s[i]) { 
		  s[i]->sample(g,step,fwd,t[i]->input,
			       w->getStateArray()[t[i]->iwaveindex],
			       t[i]->rarrindex,t[i]->iwaveindex,stream,
			       dryrun,drystr);
		}
	      }
#ifdef IWAVE_VERBOSE
	      cerr<<"backwards step"<<at<<"after sample iwdx=1 ridx=1 ucb="<<
		(w->getRDOMArray()[order]->_s)[1]._s0[48]<<"\n";
#endif
	      if (!dryrun) {
		for (int iv=fd_numsubsteps(ic)-1; iv>-1; iv--) {
		  ic.get_timestep()(w->getRDOMArray(),fwd,iv,fdm);
		  for (size_t k=pow2(order-1); k<pow2(order); k++) {
		    //		    if (w->getStateArray()[k]->printact > 5) 
		    //		      fprintf(stream,"\n*** SYNCH: iwdx = %d\n\n",k);
		    synch(w->getStateArray()[k],fwd,at,iv,ic,stream);
		  }
		}	
#ifdef IWAVE_VERBOSE
		cerr<<"backwards step"<<at<<"after time step iwdx=1 ridx=1 ucb="<<
		  (w->getRDOMArray()[order]->_s)[1]._s0[48]<<"\n";
		cerr<<"IWaveSim::run adj step "<<at;
#endif
	      }
	      else {
		drystr<<"IWaveSim::run adj step "<<at;
	      }
	      at--;
#ifdef IWAVE_VERBOSE
	      cerr<<"->"<<at<<"; ref step ="<<it<<endl;
#endif
	      if (dryrun) {
		drystr<<"->"<<at<<"; ref step ="<<it<<endl;
	      }
	    }
	    if (whatodo == ACTION::restore) {
	      int cp = r.getcheck();
	      it = cplist.at(cp);
	      for (int j=0;j<(int)pow2(order-1);j++) {
		int l = 0;
		for (int k=0;k<RDOM_MAX_NARR;k++) {
		  if (fd_isarr(k,w->getStateArray()[0]->model,ic) && fd_isdyn(k,ic)) {		  
		    if (ra_a_copy(&(((w->getRefRDOMArray())[j])->_s[k]),&(cps[cp][j][l]))) {
		      RVLException e;
		      e<<"Error: IWaveSim::run\n";
		      e<<"attempt to restore checkpoint "<<cp
		       <<" encountered error from ra_a_copy\n";
		      e<<"IWAVE index = "<<j<<" checkpoint RARR index = "
		       <<l<<" RDOM RARR index = "<<k<<"\n";
		      throw e;
		    }
		    l++;
		  }
		}
	      }
#ifdef IWAVE_VERBOSE
	      cerr<<"\nIWaveSim::run adj - restored step "<<it<<" from checkpoint "<<cp<<endl;
#endif
	      if (dryrun) {
		drystr<<"\nIWaveSim::run adj - restored step "<<it<<" from checkpoint "<<cp<<endl;
	      }
	    }
	    if (whatodo == ACTION::error) {
	      RVLException e;
	      e<<"Error: IWaveSim::run, adjoint mode\n";
	      e<<"  irregular termination of revolve\n";
	      throw e;
	    }
	  }
	  while ((whatodo != ACTION::terminate) && (whatodo != ACTION::error));

	  // final sample
	  if (dryrun) drystr<<"\n";
	  step[g.dim]=at;
	  for (size_t i=0; i<t.size(); i++) {
	    // pert sample
	    if (t[i]->iwaveindex >= (int)pow2(order-1) && s[i]) { 
	      s[i]->sample(g,step,fwd,t[i]->input,
			   w->getStateArray()[t[i]->iwaveindex],
			   t[i]->rarrindex,t[i]->iwaveindex,stream,
			   dryrun,drystr);
	    }
	  }
	  // for dry run, dump revolve output to
	  // dry run file. Else to sim output.
	  if (dryrun) drystr<<"\nRevolve report:\n"<<ostr.str()<<"\n";
	  else (fprintf(stream,"\nRevolve report:\n%s\n",ostr.str().c_str()));

	}
	// extended axis update
	// single sim case - no extended axes, only time axis
	// NOTE THIS SHOULD NOW BE TAKEN CARE OF, as FIRST=LAST
	//	if (g.gdim <= g.dim+1) more=false;
	/* why not...
	   else {
	   bool bump = false; // all dims > i are NOT full
	   for (int i=g.dim-1;i>g.dim;i--) {
	   if (step[i]<stop[i]) {
	   step[i]++;
	   bump=true;
	   }
	   else {
	   if (bump) step[i]=start[i];
	   }
	   }
	   }
	*/
	// sanity check
	if (!(next_step(g,step)) && (panelindex != panelnum-1)) {
	  RVLException e;
	  e<<"Error: IWaveSim::run\n";
	  e<<"  step array at last element but \n";
	  e<<"  panelindex = "<<panelindex<<" not =\n";
	  e<<"  panelnum   = "<<panelnum<<"\n";
	  throw e;
	}
	if (dryrun) {
	  for (int i=g.dim+1;i<g.gdim;i++) {
	    drystr<<"NEXT: step["<<i<<"]="<<step[i]<<"\n";
	  }
	}
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from IWaveSim::run\n";
      throw e;
    }
  }

  ostream & IWaveSim::write(ostream & str) const {
    str<<"IWaveSim ";
    if (fwd) str<<"forward simulator\n";
    else str<<"adjoint simulator\n";
    str<<"  derivative order = "<<order<<"\n";
    str<<"  number of dynamic arrays in each RDOM = "<<ndyn<<"\n";
    int nio = 0;
    for (size_t i=0;i<t.size();i++) {
      if (s[i]) nio++;
    }
    str<<"  number of i/o tasks = "<<nio<<"\n";
    if (!fwd) {
      str<<"  number of checkpoint states = "<<snaps<<"\n";
      str<<"  number of checkpoint RARRs  = "<<narr<<"\n";
    }
    return str;
  } 
  
  void IWaveSim::printgrid(FILE * fp) const {
    fprint_grid(fp, g);
  }

  void IWaveApply(int argc, char ** argv) {
    try {
      // set up environment
      PARARRAY * par = NULL;
      FILE * stream = NULL;
#ifdef IWAVE_VERBOSE
      cerr<<"IWaveApply -> IWaveEnvironment\n";
#endif
      IWaveEnvironment(argc, argv, 0, &par, &stream);

#ifdef IWAVE_USE_MPI
      if (retrieveGroupID() == MPI_UNDEFINED) {
#ifdef IWAVE_VERBOSE
	fprintf(stderr,"NOTE: idle rank=%d, finalize MPI, cleanup, exit\n",retrieveGlobalRank());
#endif	
	fprintf(stream,"NOTE: finalize MPI, cleanup, exit\n");
	fflush(stream);
      }
      else {
#endif
#ifdef IWAVE_VERBOSE
	fprintf(stderr,"NOTE: working rank=%d proceed with sim\n",retrieveGlobalRank());
#endif	
	int dryrun = 0;
	parse(*par,"dryrun",dryrun);
	ofstream * drystr = NULL;
	if (dryrun) {
	  stringstream dryfile;
	  dryfile <<"dryrun";
	  dryfile <<retrieveGroupID();
	  drystr = new ofstream(dryfile.str().c_str());
	}

	// read basic job params from parfile
	// these are mandatory, more or less
	int deriv=0;
	if (!parse(*par, "deriv", deriv)) {
	  RVLException e;
	  e<<"Error: IWaveApply: failed to read parameter deriv\n";
	  e<<"  (int) = order of derivative to compute\n";
	  e<<"  values: 0 = forward map;\n";
	  e<<"          1 = first derivative or adjoint first derivative\n";
	  e<<"          2 = second derivative or adjoint second derivative\n";
	  e<<"          ...\n";
	  e<<"  check parameter file\n";
	  throw e;
	}
#ifdef IWAVE_VERBOSE
	cerr<<"IWaveApply: deriv = "<<deriv<<endl;
#endif

	bool fwd = true;
	int adj = 0;
	if (!parse(*par, "adjoint", adj) && deriv>0) {
	  RVLException e;
	  e<<"Error: IWaveApply: failed to read parameter adj\n";
	  e<<"  required for application of derivative of any positive order\n";
	  e<<"  (int) = [0 = apply derivative, 1 = apply adjoint derivative]\n";
	  e<<"  check parameter file\n";
	  throw e;
	}
	if (adj == 0) fwd = true;
	else fwd = false;
#ifdef IWAVE_VERBOSE
	cerr<<"IWaveApply: fwd = "<<fwd<<endl;
#endif
	int snaps=0;
	if (!parse(*par,"nsnaps",snaps) && !fwd) {
	  RVLException e;
	  e<<"Error: IWaveApply, adjoint case: failed to read parameter nsnaps\n";
	  e<<"  (int) = number of reference field Cauchy data checkpoints\n";
	  e<<"  check parameter file\n";
	  throw e;
	}
#ifdef IWAVE_VERBOSE
	cerr<<"IWaveApply: snaps = "<<snaps<<endl;
#endif
	int printact = 0;
	parse(*par,"printact",printact);
#ifdef IWAVE_VERBOSE
	cerr<<"IWaveApply: printact = "<<printact<<endl;
#endif
	IWaveInfo ic;
#ifdef IWAVE_VERBOSE      
	cerr<<"IWaveApply -> IWaveSim constructor\n";
	cerr<<"IWaveApply: fwd="<<fwd<<" deriv="<<deriv<<endl;
#endif
	if (dryrun) {
	  IWaveSim sim(deriv,fwd,*par,stream,ic,printact,snaps,dryrun,*drystr);
#ifdef IWAVE_VERBOSE
	  cerr<<"IWaveApply -> IWaveSim::run\n";
#endif
	  sim.run();
	}
	else {
	  IWaveSim sim(deriv,fwd,*par,stream,ic,printact,snaps,dryrun,cerr);
#ifdef IWAVE_VERBOSE
	  cerr<<"IWaveApply -> IWaveSim::run\n";
#endif
	  sim.run();
	}
#ifdef IWAVE_VERBOSE
	cerr<<"IWaveApply -> exit\n";
#endif
	if (drystr) {
	  drystr->close();
	  delete drystr;
	}
#ifdef IWAVE_USE_MPI
	/* end nontriv comm branch */
      }
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      ps_delete(&par);
      fclose(stream);
    }
    catch (RVLException & e) {
      e<<"\ncalled from IWaveApply\n";
      throw e;
    }
  }

} // end namespace TSOpt

