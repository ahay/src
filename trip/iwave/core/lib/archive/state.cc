#include "state.hh"

//#define DEBUG
#undef LALA
//#define LALA
//#define VERBOSE

namespace TSOpt {

  void IWaveEnvironment(int argc, char ** argv, int ts,
			PARARRAY ** par, 
			FILE ** stream) {
    int err=0;
    initparallel_global(ts);
    if ((err=initoutstream(stream,retrieveGlobalRank(),retrieveGlobalSize()))) {
      RVLException e;
      e<<"Error: IWaveEnvironment from initoutstream, err="<<err<<"\n";
      e<<"failed to initialize output stream\n";
      throw e;
    }
    if ((err=readinput(par,*stream,argc,argv))) {
      RVLException e;
      e<<"Error: IWaveEnvironment from readinput, err="<<err<<"\n";
      e<<"failed to parse param file\n";
      throw e;
    }
    if ((err=initparallel_local(*(*par),*stream))) {
      RVLException e;
      e<<"Error: IWaveEnvironment from initparallel_local, err="<<err<<"\n";
      e<<"failed to create cartesian grid or local comm\n";
      throw e;
    }
  }

  IWaveState::IWaveState(PARARRAY  & _pars, FILE * _stream,
			 GFDM_INIT_FUN _minit)
      : stream(_stream), pars(_pars), isinit(false), tsi(iwstate.model.tsind)
       {
    _minit(&gfdm);
  }

  IWaveState::~IWaveState() { 
    //    cerr<<"IWaveState::~IWaveState, call iwave_destroy\n";
    if (isinit) iwave_destroy(&iwstate); 
    //    cerr<<"IWaveState::~IWaveState, exit\n";
  }

  void IWaveState::setTime(Time const & t) {
    try {
      TSIndex const & ts1 = dynamic_cast< TSIndex const &>(t);
      TIMESTEPINDEX & tin = (this->tsi).getCstruct();
      tin.it=ts1.getCstruct().it;
      tin.iv=ts1.getCstruct().iv;
      tin.dt=ts1.getCstruct().dt;
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: IWaveState::setTime\n";
      e<<"input Time object not TSIndex\n";
      throw e;
    }
  }

  void IWaveState::initstate() const {

    // GFD_MODEL data member initialized on construction!
    // use IWAVE::FD_MODEL intialization function to initialize IWAVE object
    // int err=iwave_construct(&iwstate,&pars,stream,gfdm.fd_model_init);
    // fish fd_model_init function out of gfd_model struct
    int err=iwave_construct(&iwstate,&pars,stream,gfdm.fd_model_init);

    if (err) {
      RVLException e;
      e<<"Error: IWaveState constructor, from iwave_construct, err="<<err<<"\n";
      throw e;
    }

    isinit=true;

    // print according to flags set in param table
    iwave_printf(&iwstate,&pars,stream);
  }
  
  IWAVE const & IWaveState::getIWAVE() const {
    if (!isinit) initstate();
    return iwstate;
  }

  IWAVE & IWaveState::getIWAVE() {
    if (!isinit) initstate();
    return iwstate;
  }

  IWaveLinState::IWaveLinState(PARARRAY & _pars, FILE * _stream,
			       //			       int (*gminit)(GFD_MODEL * mdf))
			       GFDM_INIT_FUN _minit)
    : IWaveState(_pars,_stream,_minit), 
      isinit(false),
      ltsi(linstate.model.tsind) {  }

  IWaveLinState::~IWaveLinState() {
    if (isinit) iwave_destroy(&linstate); 
  }

  void IWaveLinState::initstate() const {

    /* first initialize the parent if necessary */
    if (!IWaveState::isinit) IWaveState::initstate();

    int err=0;

    err=iwave_construct(&linstate,&pars,stream,gfdm.fd_model_init);
    //cerr<<"IWaveLinState::IWaveLinState,after iwave_construct, err ="<<err<<"\n";
    //cin>>err;
    if (err) {
      RVLException e;
      e<<"Error: IWaveLinState constructor, err="<<err<<"\n";
      throw e;
    }

    FD_MODEL * fdm = (FD_MODEL *)(iwstate.model.specs);
    FD_MODEL * lfdm = (FD_MODEL *)(linstate.model.specs);
    // copy parameters
    fdm->parcopy(lfdm->fdpars,fdm->fdpars);

    isinit=true;
  }

  IWAVE const & IWaveLinState::getIWAVE() const {
    if (!isinit) initstate();
    return linstate;
  }

  IWAVE & IWaveLinState::getIWAVE() {
    if (!isinit) initstate();
    return linstate;
  }    

  void IWaveLinState::setTime(Time const & t) {
    try {
      TSIndex const & ts1 = dynamic_cast< TSIndex const &>(t);
      TIMESTEPINDEX & tin = (this->ltsi).getCstruct();
      tin.it=ts1.getCstruct().it;
      tin.iv=ts1.getCstruct().iv;
      tin.dt=ts1.getCstruct().dt;
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: IWaveLinState::setTime\n";
      e<<"input Time object not TSIndex\n";
      throw e;
    }
  }

  void IWaveLinStep::run() {

    int printact = iw.getIWAVE().printact;
    int err=0;

    FILE *stream = iw.getStream();
    if(printact>1){
      fprintf(stream,"\n---enter IWaveLinStep::run() --- ");
      if(this->fwd) fprintf(stream,"fwd = true ---\n ");
      else fprintf(stream,"fwd = false ---\n"); 
    }

    int iv; int it;
    it = iw.getIWAVE().model.tsind.it;
    iv = iw.getIWAVE().model.tsind.iv;
    if (printact>5) {
      fprintf(stream,"\n*****************************************************\n");
      fprintf(stream,"\n*****************************************************\n");
      fprintf(stream,"\n*****************************************************\n");
      fprintf(stream,"* BEGIN IWaveLinStep::run(): it = %d iv = %d \n",it,iv);
      rd_a_print(&(iw.getIWAVE().model.ld_a),stream);
      fprintf(stream,"\n*****************************************************\n");
      fprintf(stream,"\n*****************************************************\n");
      fprintf(stream,"\n*****************************************************\n");
    }
    // choose between fwd lin and adj lin time steps
    GEN_TIMESTEP_FUN ts;
    if (fwd) ts = iw.getGFDM().tsfm;
    else ts = iw.getGFDM().tsam;

    // choose between fwd lin and adj line update indicators
    // also adj exchange before step in adjoint case, 
    // (sends and receives swapped), fwd exchange after step
    // in fwd case
    GEN_UPDATE_FUN ud;
#ifndef LALA
    if (fwd)     ud = iw.getGFDM().udfm;
    else {
   
      ud = iw.getGFDM().udam;
      // in adjoint case, zero receive arrays, as they will be used
      // as buffers to accumulate non-local updates
      if ((err=giwave_zero_recv(&(iw.getIWAVE()),
			       ud,
				iw.getStream()))) {
	RVLException e;;
	e<<"ERROR: IWaveLinStep::run from giwave_synch, err="<<err<<"\n";
	throw e;
      }
    }
#endif
#ifdef LALA
    ud=iw.getGFDM().udfm;
    if (!fwd) {
	if ((err=giwave_synch(&(iw.getIWAVE()),
			   ud,
			   fwd,
			      iw.getStream()))) {
	RVLException e;;
	e<<"ERROR: IWaveLinStep::run from giwave_synch, err="<<err<<"\n";
	throw e;
      }	
    }    
#endif

    if ((err=giwave_dmod(&(iw.getIWAVE()),
			&(iw.IWaveState::getIWAVE()),
			ts,
			//			ud,
			 iw.getStream()))) {
      RVLException e;;
      e<<"ERROR: IWaveLinStep::run from giwave_dmod, err="<<err<<"\n";
      throw e;
    }

#ifdef LALA
    if (fwd) {
#endif
	if ((err=giwave_synch(&(iw.getIWAVE()),
			   ud,
			   fwd,
			      iw.getStream()))) {
	RVLException e;;
	e<<"ERROR: IWaveLinStep::run from giwave_synch, err="<<err<<"\n";
	throw e;
      }	
#ifdef LALA
    }
#endif
    if (printact>5) {
      fprintf(stream,"\n*****************************************************\n");
      fprintf(stream,"\n*****************************************************\n");
      fprintf(stream,"\n*****************************************************\n");
      fprintf(stream,"* END IWaveLinStep::run(): it = %d iv = %d \n",it,iv);
      rd_a_print(&(iw.getIWAVE().model.ld_a),stream);
      fprintf(stream,"\n*****************************************************\n");
      fprintf(stream,"\n*****************************************************\n");
      fprintf(stream,"\n*****************************************************\n");
    }
    if(printact>1){
      fprintf(stream,"\n---exit IWaveLinStep::run()\n ");
    }

  }

  void IWaveStaticInit::run() {
    int err = 0;
    
    IWAVE & iw_iwave = iw.IWaveState::getIWAVE();
    FILE *stream = iw.getStream();
    
    int firstrec, lastrec;
    firstrec = lastrec = 0;
    /* compute firstpanel and lastpanel */
    calc_group(&firstrec, &lastrec, s.getNumRecs());
    
    // cerr<<"---- IWaveStaticInit::run, panelindex = "<<panelindex<<", firstpanel = "<<firstpanel<<", lastpanel = "<<lastpanel<<", panelnum = "<<panelnum<<" -----\n";
    //    fprintf(stream,"---- IWaveStaticInit::run -----\n");
    //    fprintf(stream,"---- recindex = %d, firstrec = %d, lastrec = %d, numrecs = %d -----\n", s.getRecIndx(),firstrec,lastrec,s.getNumRecs());
    if( (s.getRecIndx() >= firstrec) && (s.getRecIndx() <= lastrec) ) {
      //	  cerr<<"---- calling iwave_static_init -----"<<endl;
      //      fprintf(stream,"---- calling iwave_static_init -----\n");
      err = iwave_static_init(&(iw_iwave), 
			      &(iw.getPAR()), 
			      stream, 
			      s.getRecIndx(),
			      firstrec);     
    }
    else {
      //	  cerr<<"---- skip iwave_static_init call -----"<<endl;
      //      fprintf(stream,"---- skip iwave_static_init call-----\n");
    }
    fflush(stream);
    
    if (err) {
      RVLException e;
      e<<"ERROR: IWaveStaticInit::run \n";
      e<<"---- nonzero return = "<<err<<" from iwave_static_init -----\n";
      throw e;
    }
    //    fprintf(stream,"---- IWaveStaticInit::run, normal exit -----\n");
    fflush(stream);
  }

  void IWaveDynamicInit::setcheckTime(Time const & tin) {
    try {
      TSIndex const & tsin = dynamic_cast<TSIndex const &>(tin);
      tcheck=tsin;
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: IWaveDynamicInit::setcheckTime\n";
      e<<"input time not TSIndex object\n";
      throw e;
    }
  }

  void IWaveDynamicInit::setstartTime(Time const & tin) {
    try {
      TSIndex const & tsin = dynamic_cast<TSIndex const &>(tin);
      tstart=tsin;
      if (fwd) tcheck = tstart; /* default tcheck for forward sim */
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: IWaveDynamicInit::setstartTime\n";
      e<<"input time not TSIndex object\n";
      throw e;
    }
    
  }

  void IWaveDynamicInit::setfinalTime(Time const & tin) {
    try {
      TSIndex const & tsin = dynamic_cast<TSIndex const &>(tin);
      tfinal=tsin;
      if (!fwd) tcheck = tfinal; /* default tcheck for backward sim*/
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: IWaveDynamicInit::setstartTime\n";
      e<<"input time not TSIndex object\n";
      throw e;
    }
  }

  void IWaveDynamicInit::setTimerange(Time const & tstartin,Time const & tfinalin) {
    setstartTime(tstartin);
    setfinalTime(tfinalin);     
  }

  void IWaveDynamicInit::takeshot(){
    int err=0;
    if (iw.IWaveState::getTSIndex().getCstruct().iv == 0){
      err = iwave_dynamic_takeshot(&(iw.IWaveState::getIWAVE()),iw.IWaveState::getTSIndex().getCstruct().it-_tstart.it);
    }
    if (err) {
      RVLException e;
      e<<"Error: IWaveDynamicInit::takeshot\n";
      e<<"iwave_dynamic_takeshot return error \n";
      throw e;
    }
  }

  void IWaveDynamicInit::run() {
    int err = 0;
    if(fwd)
      err = giwave_dynamic_init(&(iw.IWaveState::getIWAVE()),_tcheck.it,_tcheck.it-_tstart.it);
    else {
      if (_tcheck.it >= _tfinal.it) {
	//  giwave_dynamic_init(&(iw.getIWAVE()),_tfinal.it,0);
	err = giwave_dynamic_init(&(iw.IWaveState::getIWAVE()),_tcheck.it,0); // enable to initialize dynamic fields at time after _tfinal
      }
      else 
	err = giwave_dynamic_init(&(iw.IWaveState::getIWAVE()),_tcheck.it,_tcheck.it-_tstart.it);
    }
    
    if (err) {
      RVLException e;
      e<<"ERROR: IWaveDynamicInit::run \n";
      e<<"---- nonzero return = "<<err<<" from giwave_dynamic_init -----\n";
      throw e;
    }
  }

  void IWaveLinDynamicInit::setTime(Time const & tin) { 
    try {
      TSIndex const & tsin = dynamic_cast<TSIndex const &>(tin);
      t=tsin; 
      isinit=true;
    }
    catch (bad_cast) {
      RVLException e;
      e<<"Error: attempt to assign non-TSIndex in IWaveLinDynamicInit::setTime\n";
      throw e;
    }
  }
  
  void IWaveLinDynamicInit::run() {
    if (!isinit) {
      RVLException e;
      e<<"Error: attempt to run IWaveLinDynamicInit before initialization of target time\n";
      throw e;
    }
    TIMESTEPINDEX tc = t.getCstruct();
    iwave_dynamic_init(&(iw.getIWAVE()),tc.it);
    iw.getIWAVE().model.tsind.iv = tc.iv;
  }

  void IWaveTargetRefTime::run() {
    try {
      //      TSIndex const & tmp_curt = dynamic_cast<TSIndex const &>(ts.getState().IWaveState::getTime());      
      TSIndex const & tmp_curt = dynamic_cast<TSIndex const &>(ts.getTime());
      TSIndex _curt(tmp_curt);
      //int ndim = ts.getState().getIWAVE().model.g.dim;
      //_curt.getCstruct().iv = 2 * ndim -1;
      _curt = tmp_curt;
      //      ts.getState().getGFDM().refsubstep(&(_curt.getCstruct().it), &(_curt.getCstruct().iv),&(ts.getState().IWaveState::getIWAVE().model));
      if (fwd) 
	ts.getState().getGFDM().linrefsubstep(&(_curt.getCstruct().it), &(_curt.getCstruct().iv),&(ts.getState().getIWAVE().model));
      else 
	ts.getState().getGFDM().adjrefsubstep(&(_curt.getCstruct().it), &(_curt.getCstruct().iv),&(ts.getState().getIWAVE().model));
#ifdef VERBOSE
      fprintf(ts.getState().getStream(),"IWaveTargetRefTime:: run - target time is\n");
      fprintf(ts.getState().getStream(),"-- it=%d iv=%d\n", _curt.getCstruct().it, _curt.getCstruct().iv);
      //      fprintf(stderr,"IWaveTargetRefTime:: run - target time is\n");
      //      fprintf(stderr,"-- it=%d iv=%d\n", _curt.getCstruct().it, _curt.getCstruct().iv);
#endif
      ref.setTargetTime(_curt);
#ifdef VERBOSE
      fprintf(ts.getState().getStream(),"IWaveTargetRefTime:: run - before first call\n");
      TSIndex const & tmp_dump = dynamic_cast<TSIndex const &>(ref.getState().IWaveState::getTime());
      fprintf(ts.getState().getStream(),"IWaveTargetRefTime:: run - start time is ");
      fprintf(ts.getState().getStream(),"-- it=%d iv=%d\n", tmp_dump.getCstruct().it, tmp_dump.getCstruct().iv);
#endif
      ref.run();
#ifdef VERBOSE
      fprintf(ts.getState().getStream(),"IWaveTargetRefTime:: first call to run - finish time is ");
      fprintf(ts.getState().getStream(),"-- it=%d iv=%d\n", tmp_dump.getCstruct().it, tmp_dump.getCstruct().iv);

#endif
      /* this second ref.run required because CPSim returns
	 prematurely from restore loop for some reason - should be
	 harmless in any case, but should not be required. To-Do: fix
	 CPSim so that full update occurs in one call
      */
      ref.run();
#ifdef VERBOSE
      fprintf(ts.getState().getStream(),"IWaveTargetRefTime:: second call to run - finish time is ");
      fprintf(ts.getState().getStream(),"-- it=%d iv=%d\n", tmp_dump.getCstruct().it, tmp_dump.getCstruct().iv);
#endif
    }
    catch (RVLException & e) {
      e<<"\ncalled from IWaveTargetRefTime::run\n";
      throw e;
    }
  }
  


}
