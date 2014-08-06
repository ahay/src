/** \file state.hh
 * IWAVE++ state and auxiliary algorithm classes
 * 
 * Contains the IWaveState and IWaveLinState classes, 
 * and some auxiliary algorithm classes to initialize states, 
 * run time-stepping step, store dynamic fields of states,
 * and update state time, etc.
 *
 * \author William W. Symes and Dong Sun     
 */
/*================================================================================*/
#ifndef __IWAVEPP_STATE
#define __IWAVEPP_STATE

/*--------------------------------------------------------------------------------*/
#include "except.hh"
#include "tsind.hh"
#include "sim.hh"
#include "samp.hh"
#include "seamx_headers.hh"
#include "giwave.h"

/*--------------------------------------------------------------------------------*/
namespace TSOpt {

  using RVL::RVLException;

  /** Function to initialize parallel environment, output stream, and param table */
  void IWaveEnvironment(int argc, char ** argv, int ts, PARARRAY ** par, FILE ** stream);

  /** \brief IWaveState class for reference state
   *
   * This class holds all static (material parameters) fields, 
   * all dynamic wavefields, and time, and provides methods to 
   * get access to those information appropriately. 
   * 
   */

  class IWaveState {

  private:
    
    IWaveState();
    IWaveState(IWaveState const &);    

  protected:

    FILE * stream;               /**< output stream */
    PARARRAY & pars;             /**< parameter array */
    mutable IWAVE iwstate;       /**< structure for all fields and parallel info */
    mutable bool isinit;         /**< initilization flag */
    TSIndex tsi;                 /**< structure for time info - view of iwstate 
				    TIMESTEPINDEX struct member */
    void initstate() const;      /**< to initialize IWAVE data member */

  protected:

    GFD_MODEL gfdm;
    
  public:
    
    IWaveState(PARARRAY & _pars, FILE * _stream,
	       GFDM_INIT_FUN _minit);
    virtual ~IWaveState();
    
    virtual void setTime(Time const & t);
    virtual Time & getTime() { return tsi; }
    virtual Time const & getTime() const {return tsi; }

    virtual TSIndex const & getTSIndex() const { return tsi; }
    virtual TSIndex & getTSIndex() { return tsi; }

    /** these functions can be overridden by child IWaveLinState
	to return the tangent vector instead of the base */
    virtual IWAVE const & getIWAVE() const;
    virtual IWAVE & getIWAVE();
    
    /** const access to GFDM needed at C interface arguments */
    GFD_MODEL const & getGFDM() const { return gfdm; }
    //    GFD_MODEL const & getGFDM() const { return gfdm; }
    //    FDM_INIT_FUN getFDMINIT() const { return minit; }
    //    GFDM_INIT_FUN getGFDMINIT() const { return minit; }

    PARARRAY const & getPAR() const { return pars; }
    PARARRAY & getPAR() { return pars; }
    FILE * getStream() const { return stream; }
    
    virtual ostream & id(ostream & str) const {return str<<"IWaveState"; }
    void cid(FILE * str) const { fprintf(str,"IWaveState"); }
    void printpars() const { ps_printall(pars,stream); }
    void printf() const { iwave_printf(&iwstate,&pars,stream); }

    virtual ostream & write(ostream & str) const {
      str<<"IWaveState: it = "<<IWaveState::getTSIndex().getCstruct().it<<
	", iv ="<<IWaveState::getTSIndex().getCstruct().iv<<
	", dt ="<<IWaveState::getTSIndex().getCstruct().dt<<"\n";
      
      return str;
    }
    
  };

  class IWaveLinState: public IWaveState {

  private:
    
    mutable IWAVE linstate;
    mutable bool isinit;
    TSIndex ltsi;
    void initstate() const;
    IWaveLinState();
    IWaveLinState(IWaveLinState const &);
    
  public:
    
    IWaveLinState(PARARRAY & _pars, FILE * _stream,
		  //		  int (*gminit)(GFD_MODEL * mdl));
		  GFDM_INIT_FUN _minit);
    virtual ~IWaveLinState();

    /** overrides */
    void setTime(Time const & t);
    Time & getTime() { return ltsi; }
    Time const & getTime() const {return ltsi; }

    TSIndex const & getTSIndex() const { return ltsi; }
    TSIndex & getTSIndex() { return ltsi; }

    IWAVE const & getIWAVE() const;
    IWAVE & getIWAVE();
    
    ostream & id(ostream & str) const {return str<<"IWaveLinState"; }
    void cid(FILE * str) const { fprintf(str,"IWaveLinState"); }
    
    ostream & write(ostream & str) const {
      str<<"IWaveLinState: it = "<<IWaveLinState::getTSIndex().getCstruct().it<<
	", iv ="<<IWaveLinState::getTSIndex().getCstruct().iv<<
	", dt ="<<IWaveLinState::getTSIndex().getCstruct().dt<<"\n";
      
      return str;
    }

  };

  class IWaveStep: public TimeStep<IWaveState> {
  private:
    IWaveState & iw;
    IWaveStep();
    IWaveStep(IWaveStep const &);
  public:
    IWaveStep(IWaveState & _iw): iw(_iw) {}
    ~IWaveStep() {}
    IWaveState & getState() { return iw; }
    IWaveState const & getState() const { return iw; }    
    Time const & getTime() const { return iw.IWaveState::getTime(); }
    /// LOOK OUT!!!
    Time const & getNextTime() const { 
      RVLException e;
      e<<"Error: IWaveStep::getNextTime not implemented\n";
      throw e;
    }

    void setTime(Time const & t) { iw.IWaveState::setTime(t); }

    void run() {

      /*      
      fprintf(iw.getStream(),"\n**********enter IWaveStep::run for ");
      iw.cid(iw.getStream());
      fprintf(iw.getStream(),"\n");
      */
      
      int err=0;
      err=iwave_run(&(iw.IWaveState::getIWAVE()),iw.IWaveState::getStream());

      if (err) {
	RVLException e;
	e<<"Error: IWaveStep from giwave_run, err="<<err<<"\n";
	throw e;
      }

      /*      
      fprintf(iw.getStream(),"\n********** exit IWaveStep::run for ");
      iw.cid(iw.getStream());
      fprintf(iw.getStream(),"\n");
      */
    }

    ostream & write(ostream & str) const {
      str<<"IWaveStep\n";
      return str;
    }
  };

  /** linearization with respect to model, or "Born source", step -
      forward and adjoint */
  class IWaveLinStep: public TimeStep<IWaveLinState> {
    
  private:
    
    IWaveLinState & iw;
    const bool fwd;
  
    IWaveLinStep();
    IWaveLinStep(IWaveLinStep const &);

  public: 
    
    IWaveLinStep(IWaveLinState & _iw, const bool _fwd = true)
      : iw(_iw), fwd(_fwd) {}

    IWaveLinState & getState() { return iw; }
    IWaveLinState const & getState() const { return iw; }    
    Time const & getTime() const { return iw.IWaveLinState::getTime(); }
    /// LOOK OUT!!!
    Time const & getNextTime() const { 
      RVLException e;
      e<<"Error: IWaveStep::getNextTime not implemented\n";
      throw e;
    }
    void setTime(Time const & t) { iw.IWaveLinState::setTime(t); }
    void run();
 
    ostream & write(ostream & str) const {
      str<<"IWaveLinStep\n";
      return str;
    }
  };

  /** Initialization alg for static data of state (coefficient arrays) */
  class IWaveStaticInit: public Algorithm {
  private:
    IWaveState & iw;                   /* nonconst reference to external state */
    Sampler const & s;                 /* const reference to sampler coord object */
    IWaveStaticInit();
    IWaveStaticInit(IWaveStaticInit const &);
  public:
    IWaveStaticInit(IWaveState & _iw, Sampler const & _s)
      : iw(_iw), s(_s) {}
    ~IWaveStaticInit() {}
    void run();
  };

  /** Initialization alg for dynamic data of state */
  class IWaveDynamicInit: public Algorithm {
  private:
    IWaveState & iw;        /* nonconst reference to external state */
    const bool fwd;       /* true, forward; false, backward */
    mutable TIMESTEPINDEX _tstart, _tcheck, _tfinal;    
    mutable TSIndex tstart, tcheck, tfinal;
    IWaveDynamicInit();
    IWaveDynamicInit(IWaveDynamicInit const &);
  public:
    IWaveDynamicInit(IWaveState & _iw, const bool _fwd = true)
      : iw(_iw), fwd(_fwd), tstart(_tstart), tcheck(_tcheck), tfinal(_tfinal) {}
    ~IWaveDynamicInit() {}
    
    void setcheckTime(Time const & tin);
    void setstartTime(Time const & tin);
    void setfinalTime(Time const & tin);   
    void setTimerange(Time const & tstartin,Time const & tfinalin);
    void takeshot();
    void run();
  };

  /** Initialization alg for dynamic data of linearized state */
  class IWaveLinDynamicInit: public Algorithm {
  private:
    IWaveLinState & iw;        /* nonconst reference to external state */
    mutable TSIndex t;
    bool isinit;
    IWaveLinDynamicInit();
    IWaveLinDynamicInit(IWaveLinDynamicInit const &);
  public:
    IWaveLinDynamicInit(IWaveLinState & _iw )
      : iw(_iw), isinit(false) {}
    ~IWaveLinDynamicInit() {}

    void setTime(Time const & tin);
      
    void run();
  };


  /** written as template to accomodate both state and lin state */
  template<typename State>
  class IWaveNextTime: public Algorithm {
  private:
    State & iw;
    const bool fwd;  /* true, forward; false, backward */
    int niv;         /* number of substeps */
    IWaveNextTime();
    IWaveNextTime(IWaveNextTime<State> const &);
  public:
    IWaveNextTime(State & _iw, const bool _fwd = true): iw(_iw), fwd(_fwd) {
      FD_MODEL * fdm = (FD_MODEL *)(iw.getIWAVE().model.specs);
      niv = fdm->numsubsteps();
    }
    ~IWaveNextTime() {}
    void run() {
      TIMESTEPINDEX & ts = iw.State::getTSIndex().getCstruct();
      //      cerr<<"--- iwavenexttime enter with state ";
      //      iw.State::getTSIndex().write(cerr);
      //      fprintf(iw.State::getStream(),"--- iwavenexttime, enter with it=%d iv=%d\n",ts.it,ts.iv);
      //cerr<<"--- iwavenexttime, enter with it "<< ts.it <<", iv "<<ts.iv<<endl;
      if(fwd) {
	if (ts.iv < niv-1 && ts.iv >= 0)
	  ts.iv++;
	else if (ts.iv == niv-1){
	  ts.iv = 0;
	  ts.it++;
	}
	else {
	  cerr<<"Error: IWaveNextTime::run, internal action out of range \n";
	  RVLException e;
	  e<<"Error: IWaveNextTime::run, internal action out of range \n";
	  throw e;
	}
      }
      else {
	if (ts.iv <= niv-1 && ts.iv > 0)
	  ts.iv--;
	else if (ts.iv == 0){
	  ts.iv = niv - 1;
	  ts.it--;
	}
	else {
	  cerr<<"Error: IWaveNextTime::run, internal action out of range \n";
	  RVLException e;
	  e<<"Error: IWaveNextTime::run, internal action out of range \n";
	  throw e;
	}
      }
      iw.State::getIWAVE().model.tsind.it = ts.it;
      iw.State::getIWAVE().model.tsind.iv = ts.iv;
      //      cerr<<"--- iwavenexttime exit with state ";
      //      iw.State::getTSIndex().write(cerr);
      //      fprintf(iw.State::getStream(),"--- iwavenexttime, exit with it=%d iv=%d\n",ts.it,ts.iv);
    }
  };

  /** runs member Sim (reference simulator) until it reaches the
      target refference time level of the member TimeStep. */
  class IWaveTargetRefTime: public Algorithm {
  private:
    Sim<IWaveState> & ref;
    IWaveLinStep const & ts;
    bool fwd; 
    
    IWaveTargetRefTime();
  public:
    IWaveTargetRefTime(Sim<IWaveState> & _ref, 
		       IWaveLinStep const & _ts,
		       bool _fwd)
      : ref(_ref), ts(_ts), fwd(_fwd) {}

    IWaveTargetRefTime(IWaveTargetRefTime const & t)
      : ref(t.ref), ts(t.ts), fwd(t.fwd) {}
    ~IWaveTargetRefTime() {}

     void run();
  };

}

#endif
