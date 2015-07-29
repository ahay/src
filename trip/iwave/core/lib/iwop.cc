#include "iwop.hh"

namespace TSOpt {

  IWaveSpace::IWaveSpace(PARARRAY const & par, 
			 IWaveInfo const & ic,
			 bool input,
			 ostream & outfile): _s(0), _keys(0) {
    try {

      for (int j=0;j<ic.get_num_iokeys();j++) {

	if ((input == ic.get_iwave_iokeys()[j].input) &&
	    ic.get_iwave_iokeys()[j].active) {
	    
	  string hdr = "";
	  string key = ic.get_iwave_iokeys()[j].keyword;
	  if (!parse(par,key,hdr)) {
	    RVLException e;
	    e<<"Error: IWaveSpace constructor\n";
	    e<<"  key "<<ic.get_iwave_iokeys()[j].keyword<<" does not match\n";
	    e<<"  any in param array\n";
	    throw e;
	  }

	  // extract suffix - this is the obvious "case"
	  // implementation, but unsatisfactory in that adding a new
	  // space type requires hacking the code. Something more
	  // like a policy would be better, but not obvious how to
	  // carry that off
	  size_t pos = hdr.find_last_of(".");
	  if (pos==std::string::npos || pos >= hdr.size()-1) {
	    RVLException e;
	    e<<"Error: IWaveSpace constructor - filename "<<hdr<<" has no suffix\n";
	    throw e;
	  }
	  size_t net = hdr.size()-pos-1;
	  string suffix = hdr.substr(pos+1,net);	    
	    
	  // case list
	  Space<ireal> * sp = NULL;
	    
	  if (suffix == "rsf") {
#ifdef IWAVE_USE_MPI
	    sp = new MPIGridSpace(hdr,key,true,retrieveGlobalComm(),outfile);
#else 
	    sp = new GridSpace(hdr,key,true,outfile);	      
#endif
	  }
	  else if (suffix == "su") {
#ifdef IWAVE_USE_MPI
	    sp = new MPISEGYSpace(hdr,key,retrieveGlobalComm(),outfile);
#else 
	    sp = new SEGYSpace(hdr,key,outfile);	      
#endif	      
	  }
	  if (sp) {
	    _s.push_back(sp);
	    _keys.push_back(key);
	  }	    
	}
      }
      if (_s.size()==0) {
	RVLException e;
	e<<"Error: IWaveSpace constructor\n";
	e<<"  failed to locate any spaces\n";
	throw e;
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from IWaveSpace constructor\n";
      throw e;
    }
  }

  IWaveSpace::IWaveSpace(IWaveSpace const & sp)
    : _s(sp._s.size()), _keys(sp._keys) {
    try {
      for (size_t i=0;i<sp._s.size(); i++) {
	// case list, again - icky poo
#ifdef IWAVE_USE_MPI
	MPIGridSpace const * gsp = NULL;
	MPISEGYSpace const * tsp = NULL;
	if (gsp=dynamic_cast<MPIGridSpace const *>((sp._s)[i])) {
	  _s[i]=new MPIGridSpace(*gsp);
	}
	else if (tsp=dynamic_cast<MPISEGYSpace const *>((sp._s)[i])) {
	  _s[i]=new MPISEGYSpace(*tsp);	    
	}
	else {
	  RVLException e;
	  e<<"Error: IWaveSpace copy constructor\n";
	  e<<"  element "<<i<<" of prototype is neither MPIGridSpace nor MPISEGYSpace\n";
	  throw e;
	}
#else
	GridSpace const * gsp = NULL;
	SEGYSpace const * tsp = NULL;
	if ((gsp=dynamic_cast<GridSpace const *>((sp._s)[i]))) {
	  _s[i]=new GridSpace(*gsp);
	}
	else if ((tsp=dynamic_cast<SEGYSpace const *>((sp._s)[i]))) {
	  _s[i]=new SEGYSpace(*tsp);	    
	}
	else {
	  RVLException e;
	  e<<"Error: IWaveSpace copy constructor\n";
	  e<<"  element "<<i<<" of prototype is neither GridSpace nor SEGYSpace\n";
	  throw e;
	}
#endif
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from IWaveSpace copy constructor\n";
      throw e;
    }
  }
    
  IWaveSpace::~IWaveSpace() {
    for (size_t i=0;i<_s.size();i++) {
      delete _s[i];
    }
  }

  /** implements virtual DataContainer constructor via
      StdProductDataContainer class. */
  DataContainer * IWaveSpace::buildDataContainer() const {
    StdProductDataContainer * d = 
      new StdProductDataContainer();
    for (size_t i=0; i<getSize(); i++) {
      SpaceDCF<ireal> f(*(_s[i]));
      d->push(f);
    }
    return d;
  }

  size_t IWaveSpace::getSize() const { return _s.size(); }
  Space<ireal> const & IWaveSpace::operator[](size_t i) const { return *(_s[i]); }
  std::vector<std::string> IWaveSpace::getKeys() const { return _keys; }

  void IWaveOp::param_set(RVL::Vector<ireal> const & x, 
			  PARARRAY & pars, 
			  IWaveSpace const & sp,
			  std::string const & suf,
			  FILE * stream) const {
    try {
      Components<ireal> cx(x);
      if (cx.getSize() != sp.getSize()) {
	RVLException e;
	e<<"Error: IWaveOp::param_set\n";
	e<<"  vector, key list have different numbers of components\n";
	e<<"  key list:\n";
	for (size_t i=0;i<sp.getKeys().size();i++) e<<"    keys[i]="<<sp.getKeys()[i]<<"\n";
	if (suf.size()>0) e<<"  suffix = "<<suf<<"\n";
	throw e;
      }
      for (size_t i=0;i<sp.getKeys().size();i++) {
	AssignParams ap(pars,sp.getKeys()[i]+suf,stream);
	cx[i].eval(ap);
      }
    }
    catch (RVLException & e) {
      e<<"\ncalled from IWaveOp::param_set\n";
      throw e;
    }
  }

  void IWaveOp::apply(const Vector<ireal> & x, 
		      Vector<ireal> & y) const {

    try {

      SpaceTest(this->getDomain(),x,"TSOpt::IWaveOp::apply (dom)");
      SpaceTest(this->getRange(),y,"TSOpt::IWaveOp::apply (rng)");

      PARARRAY * locpars = ps_new();
      ps_copy(&locpars,*pars);
      std::string nix = "";
      param_set(x,*locpars,dom,nix,stream);
      param_set(y,*locpars,rng,nix,stream);

      // Change 11.05.10 - print par file used in sim
      if (dump_pars) {
	fprintf(stream,"PARAM ARRAY CREATED IN IWOP::APPLY\n");
	ps_printall(*locpars,stream);
	fflush(stream);
      }

      /* zero output */
      y.zero();

      if (dump_steps) {
	fprintf(stream,"IWOP::APPLY -> STEP\n");
	fflush(stream);
      }

      int deriv=0;
      bool fwd=true;
      int snaps=0;
      int printact=0;
      parse(*locpars,"printact",printact);
      IWaveSim sim(deriv,fwd,*locpars,stream,ic,printact,snaps,dryrun,drystr,announce);
      sim.run();

      ps_delete(&locpars);

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
      
  void IWaveOp::applyDeriv(const Vector<ireal> & x, 
			   const Vector<ireal> & dx,
			   Vector<ireal> & dy) const {

    try {

      // sanity test arguments
      SpaceTest(this->getDomain(),x,"TSOpt::IWaveOp::applyDeriv (dom, ref)");
      SpaceTest(this->getDomain(),dx,"TSOpt::IWaveOp::applyDeriv (dom, pert)");
      SpaceTest(this->getRange(),dy,"TSOpt::IWaveOp::applyDeriv (rng)");

      PARARRAY * locpars = ps_new();
      ps_copy(&locpars,*pars);
      std::string nix = "";
      std::string dsuf = "_d1";
      param_set(x,*locpars,dom,nix,stream);
      param_set(dx,*locpars,dom,dsuf,stream);
      param_set(dy,*locpars,rng,nix,stream);

      if (dump_pars) {
	fprintf(stream,"PARAM ARRAY CREATED IN IWOP::APPLYDERIV\n");
	ps_printall(*locpars,stream);
	fflush(stream);
      }

      if (dump_steps) {
	fprintf(stream,"IWOP::APPLYDER -> STEP\n");
	fflush(stream);
      }

      int deriv=1;
      bool fwd=true;
      int snaps=0;
      int printact=0;
      parse(*locpars,"printact",printact);

      IWaveSim sim(deriv,fwd,*locpars,stream,ic,printact,snaps,dryrun,drystr,announce);
      sim.run();

      ps_delete(&locpars);

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
      
  void IWaveOp::applyAdjDeriv(const Vector<ireal> & x, 
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
	
      PARARRAY * locpars = ps_new();
      ps_copy(&locpars,*pars);
      std::string nix = "";
      std::string dsuf = "_b1";
      param_set(x,*locpars,dom,nix,stream);
      param_set(dx,*locpars,dom,dsuf,stream);
      param_set(dy,*locpars,rng,nix,stream);

      if (dump_pars) {
	fprintf(stream,"PARAM ARRAY CREATED IN IWOP::APPLYADJDERIV\n");
	ps_printall(*locpars,stream);
	fflush(stream);
      }

      int deriv=1;
      bool fwd=false;
	
      int snaps=DEFAULT_SNAPS;
      parse(*locpars,"nsnaps",snaps);
      int printact=0;
      parse(*locpars,"printact",printact);

      if (dump_steps) {
	fprintf(stream,"IWOP::APPLYADJDER -> BUILD IWAVESIM\n");
	fflush(stream);
      }      

      IWaveSim sim(deriv,fwd,*locpars,stream,ic,printact,snaps,dryrun,drystr,announce);
      if (dump_steps) {
	fprintf(stream,"IWOP::APPLYADJDER -> RUN IWAVESIM\n");
	fflush(stream);
      }      
      sim.run();

      ps_delete(&locpars);

      if (dump_steps) {
	fprintf(stream,"IWOP::APPLYADJ -> BARRIER\n");
	fflush(stream);	
      }

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
      
  void IWaveOp::applyDeriv2(const Vector<ireal> & x, 
			    const Vector<ireal> & dx0,
			    const Vector<ireal> & dx1,
			    Vector<ireal> & dy) const {

    try {

      // sanity test arguments
      SpaceTest(this->getDomain(),x,"TSOpt::IWaveOp::applyDeriv2 (dom, ref)");
      SpaceTest(this->getDomain(),dx0,"TSOpt::IWaveOp::applyDeriv2 (dom, pert1)");
      SpaceTest(this->getDomain(),dx1,"TSOpt::IWaveOp::applyDeriv2 (dom, pert2)");
      SpaceTest(this->getRange(),dy,"TSOpt::IWaveOp::applyDeriv2 (rng)");

      PARARRAY * locpars = ps_new();
      ps_copy(&locpars,*pars);
      std::string nix = "";
      param_set(x,*locpars,dom,nix,stream);
      std::string dsuf = "_d1";
      param_set(dx0,*locpars,dom,dsuf,stream);
      dsuf = "_d2";
      param_set(dx1,*locpars,dom,dsuf,stream);
      param_set(dy,*locpars,rng,nix,stream);

      if (dump_pars) {
	fprintf(stream,"PARAM ARRAY CREATED IN IWOP::APPLYDERIV\n");
	ps_printall(*locpars,stream);
	fflush(stream);
      }

      if (dump_steps) {
	fprintf(stream,"IWOP::APPLYDER -> STEP\n");
	fflush(stream);
      }

      int deriv=2;
      bool fwd=true;
      int snaps=0;
      int printact=0;
      parse(*locpars,"printact",printact);

      IWaveSim sim(deriv,fwd,*locpars,stream,ic,printact,snaps,dryrun,drystr,announce);
      sim.run();

      ps_delete(&locpars);

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
      
  void IWaveOp::applyAdjDeriv2(const Vector<ireal> & x, 
			       const Vector<ireal> & dx0,
			       const Vector<ireal> & dy,
			       Vector<ireal> & dx1) const {

    try{

      if (dump_steps) {
	fprintf(stream,"IWOP::APPLYADJ2 -> START\n");
	fflush(stream);
      }	

      SpaceTest(this->getDomain(),x,"TSOpt::IWaveOp::applyAdjDeriv2 (dom, ref)");
      SpaceTest(this->getDomain(),dx0,"TSOpt::IWaveOp::applyAdjDeriv2 (dom, pert1)");
      SpaceTest(this->getRange(),dy,"TSOpt::IWaveOp::applyAdjDeriv2 (rng)");
      SpaceTest(this->getDomain(),dx1,"TSOpt::IWaveOp::applyAdjDeriv2 (dom, pert2)");
	
      std::string nix = "";
      PARARRAY * locpars = ps_new();
      ps_copy(&locpars,*pars);
      param_set(x,*locpars,dom,nix,stream);
      std::string dsuf = "_d1";
      param_set(dx0,*locpars,dom,dsuf,stream);
      dsuf = "_b2";
      param_set(dx1,*locpars,dom,dsuf,stream);
      param_set(dy,*locpars,rng,nix,stream);

      if (dump_pars) {
	fprintf(stream,"PARAM ARRAY CREATED IN IWOP::APPLYADJDERIV2\n");
	ps_printall(*locpars,stream);
	fflush(stream);
      }

      if (dump_steps) {
	fprintf(stream,"IWOP::APPLYADJDER2 -> STEP\n");
	fflush(stream);
      }

      int deriv=2;
      bool fwd=false;
	
      int snaps=DEFAULT_SNAPS;
      parse(*locpars,"nsnaps",snaps);
      int printact=0;
      parse(*locpars,"printact",printact);

      IWaveSim sim(deriv,fwd,*locpars,stream,ic,printact,snaps,dryrun,drystr,announce);
      sim.run();

      ps_delete(&locpars);

#ifdef IWAVE_USE_MPI
      MPI_Barrier(retrieveGlobalComm());
#endif
      if (dump_steps) {
	fprintf(stream,"IWOP::APPLYADJDERIV2 -> EXIT\n");
	fflush(stream);	
      }
    }
    catch (RVLException & e) {
      e<<"\n called from TSOpt::IWaveOp::applyAdjDeriv\n";
      throw e;
    }

  }

  Operator<ireal> * IWaveOp::clone() const { 
    return new IWaveOp(*this);
  }
    
  IWaveOp::IWaveOp(PARARRAY _pars, FILE * _stream,
		   bool _dryrun,
		   ostream & _drystr,
		   ostream & _announce)
    : dom(_pars,ic,true,cerr), rng(_pars,ic,false,cerr), 
      stream(_stream),
      pars(NULL),
      dump_steps(0),
      dump_pars(0),
      dump_term(0),
      dryrun(_dryrun),
      drystr(_drystr),
      announce(_announce) {
      
    int err=0;
      
    /* copy input par array to data member */
    if ((err=ps_copy(&pars,_pars))) {
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
    
  IWaveOp::IWaveOp(IWaveOp const & x) 
    : dom(x.dom), rng(x.rng), 
      stream(x.stream), 
      pars(NULL),
      dump_steps(x.dump_steps),
      dump_pars(x.dump_pars),
      dump_term(x.dump_term),
      dryrun(x.dryrun),
      drystr(x.drystr),
      announce(x.announce) {
    // cerr<<"literal copy, so pararray is copied"<<endl;
    int err=0;
    if ((err=ps_copy(&pars,*(x.pars)))) {
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

  IWaveOp::~IWaveOp() { 
    ps_delete(&pars);
    /* cerr<<"good bye!\n"; */ }
      
  // added 23.06.10 to facilitate using source as variable
  // without admitting that it's part of domain
  PARARRAY & IWaveOp::getPar() { return *pars; }
  PARARRAY const & IWaveOp::getPar() const { return *pars; }

  ostream & IWaveOp::write(ostream & str) const {
    str<<"IWOP object\n";
    return str;
  }
}


