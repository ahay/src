/*************************************************************************

Copyright Rice University, 2004, 2005, 2006, 2007, 2008
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/
#ifndef __RVL_SIM
#define __RVL_SIM

#include "alg.hh"
#include "tterm.hh"

#include <stack>
#include "dtterm.hh"
#include "cdt.hh"
#include "stackBase.hh"
#include "rt.hh"

namespace TSOpt {
	
  using RVL::Writeable;
  using RVL::RVLException;
  using RVLAlg::Algorithm;
  using RVLAlg::ListAlg;
  using RVLAlg::LoopAlg;
  using RVLAlg::StateAlg;
  using RVLAlg::Terminator;

  /** Fundamental class for time step methods. Differs from StateAlg
      in (i) exposing simulation time, by delegation to TimeState
      managed by StateAlg, and (ii) providing read-only access to
      subsequent simulation time. This last method is virtual - its
      implementation is an active constraint on subclasses. Children
      expected to take the form of wrappers around less abstract code,
      perhaps collections of functions and structs.
  */

  template<typename State>
  class TimeStep: public StateAlg<State>, public Writeable {
  public:
    virtual ~TimeStep() {}
    /** convenience interface for control of simulation time */
    virtual void setTime(Time const & t) { 
      //      cout << "In TimeStep::setTime(t)" << endl;
      //      t.write(cout);
      (this->getState()).setTime(t); 
    }

    /** convenience interface to expose current simulation time */
    virtual Time const & getTime() const { return (this->getState()).getTime(); }

    /** virtual interface to expose the next simulation time */
    virtual Time const & getNextTime() const = 0;

  };

  /** basic TimeStep construction - intended to standardize
      combination of state class, possibly owning the guts of its own
      "run" method to implement time step, and the TimeStep
      interface. Note that run and write methods are still abstract in
      this class.

      The template parameter State is very general. The only real 
      constraint is that a copy method must be defined, with signature

      void State::copy(State const &),

      with the semantics of a deep copy from the argument to the State
      on which the method is called.

  */

  template<typename State>
  class StdTimeStep: public TimeStep<State> {
  protected:
   
    State mystate;    

  public:
    StdTimeStep() { }
    StdTimeStep(StdTimeStep<State> const & s): mystate(s.mystate) { }
    virtual ~StdTimeStep() {}
    
    // void setState(State const & x) { mystate.copy(x); }
    
    State & getState() {return mystate;  }
    State const & getState() const { return mystate;  }

      
  };


  /** Peeping TimeStep -- all derived classes must implement
      a peekTargetTime(Time) method that allows the timestep class
      to access the terminator's target time. Useful for adaptive
      time-stepping
  */
  template<typename State>
  class PeekTimeStep: public StdTimeStep<State> {
  protected:
   

  public:
    PeekTimeStep() { }
    PeekTimeStep(PeekTimeStep<State> const & s) { }
    virtual ~PeekTimeStep() {}
    
    virtual void peekTT(const Time & t) = 0;

    virtual double getDt() const = 0;
    virtual void setDt(double dt) = 0;
      
  };


  /** Runs main TimeStep with or without pre- or post-operations,
      represented by additional Algs. All Time-related actions
      delegated to main TimeStep alg. Fully implemented concrete
      class, not intended for further derivation. 
  */
  template<typename State>
  class TimeStepList: public TimeStep<State> {
  private:
    TimeStep<State> & main;
    Algorithm & pre;
    Algorithm & post;

    bool preflag;
    bool postflag;
    TimeStepList();
  public:
    /** Simplest constructor - pass-through to argument TimeStep */
    TimeStepList(TimeStep<State> & a)
      : main(a), pre(a), post(a), preflag(false), postflag(false) {}

    /** TimeStep argument as main Alg, with another Alg as pre- or
	post-operation, choice determined by flag: true for preop,
	false for postop. */
    TimeStepList(TimeStep<State> & a, Algorithm & b, bool pre)
      : main(a), pre(b), post(b), preflag(pre), postflag(!pre) {}

    /** TimeStep (first arg) with pre- and post-ops */
    TimeStepList(TimeStep<State> & a, Algorithm & b, Algorithm & c) 
      : main(a), pre(b), post(c), preflag(true), postflag(true) {}

    /** copy constructor */
    TimeStepList(TimeStepList<State> const & l)
      : main(l.main), pre(l.pre), post(l.post), preflag(l.preflag), postflag(l.postflag) {}
    ~TimeStepList(){}

    
    // void setState(State const & s) { main.setState(s); }
    
    State & getState() { return main.getState(); }
    State const & getState() const { return main.getState(); }
    void setTime(Time const & t) { return main.setTime(t); }
    

    Time const & getTime() const { return main.getTime(); }
    Time const & getNextTime() const { return main.getNextTime(); }
    void run() {
      if      (( preflag) && ( postflag)) {
	ListAlg inner(main,post);
	ListAlg outer(pre,inner);
	outer.run();
      }
      
      else if ((!preflag) && ( postflag)) {
	ListAlg outer(main,post);
	outer.run();
      }
      
      else if (( preflag) && (!postflag)) {
	ListAlg outer(pre,main);
	outer.run();
      }
      
      else {
	main.run();
      }

    }

    ostream & write(ostream & str) const {
      str<<"---------------- beg TimeStepList::write ----------------\n";
      str<<"TimeStepList. preop="<<preflag<<" postop="<<postflag<<"\n";
      str<<"TimeStep:\n";
      main.write(str);
      str<<"---------------- end TimeStepList::write ----------------\n";
      return str;
    }

  };

  /** Interface for timestepping simulator built out of a TimeStep and
      a TimeTerm, with most services implemented by delegation but
      run() left pure virtual. TimeStep data member intended to
      provide simulation step. TimeTerm encapsulates arbitrary
      time-dep. actions. The structure added by this class is simply
      indirection: it is a StateAlg that delegates some but not all of
      its StateAlg attributes to a member StateAlg.
  */

  template<typename State> 
  class Sim: public StateAlg<State>, public Writeable {

  protected:

    TimeStep<State> & step;
    TimeTerm & term;
    Algorithm &initstep;

    Sim();

  public: 

    Sim(TimeStep<State> & _step,
	TimeTerm & _term,
	Algorithm &_initstep)
      : step(_step), term(_term), initstep(_initstep) {}
    
    Sim(Sim<State> const & s)
      : step(s.step), term(s.term), initstep(s.initstep) {}

    virtual ~Sim() {}

    /** inherited from StateAlg, implemented by delegation to State
	member of step. */
    
    // void setState(State const & st) { step.setState(st); }
    
    /** inherited from StateAlg, implemented by delegation to State
	member of step. */ 
    State & getState() { return step.getState(); }

    /** inherited from StateAlg, implemented by delegation to State
	member of step. */ 
    State const & getState() const { return step.getState(); }

    /** convenience function: set current simulation time */
    virtual void setTime(Time const & t) { step.setTime(t); }

    /** convenience function: expose current simulation time */
    virtual Time const & getTime() const { return step.getTime(); }

    /** set target simulation time */
    virtual void setTargetTime(Time const & t) { term.setTargetTime(t); }
    
    /** get target simulation time */
    virtual Time const & getTargetTime() const { return term.getTargetTime(); }
    
    /** get initstep. Allows to create simulators that can run initialization 
        step separately (i.e. only once) from time step. Must be non-const, since
        normally modifies members of state */
    Algorithm &getInitStep() { return initstep; }

    /**  virtual function that "resets" sim object. useful for inversions. 
	 natural implementation is to run initstep, but could vary */
    virtual void resetSim() { initstep.run(); }

    /* run still virtual! */

    ostream & write(ostream & str) const {
      str<<"---------------- beg Sim::write ----------------\n";
      str<<"TimeStep:\n";
      step.write(str);
      str<<"TimeTerm:\n";
      term.write(str);
      str<<"---------------- end Sim::write ----------------\n";    
      return str;
    }
  };

  /** Argument struct for StdSim constructor - enables use of
      standard policy constructor 
  */
  template<typename State>
  class StdSimData {
  private:
    TimeStep<State> & step;
    TimeTerm & term;
    Algorithm & init;

    /** extra data used during CPSim construction */
    stackBase<State> * ldclist;
    int _steps; 
    int _snaps; 
	
    StdSimData();
  public:
    StdSimData(StdSimData<State> const & s)
      : step(s.step), term(s.term), init(s.init), ldclist(s.ldclist), 
	_steps(s._steps), _snaps(s._snaps){}
  
    /** data argument for building CPSim */
    StdSimData(TimeStep<State> & _step,
	       TimeTerm & _term,
	       Algorithm & _init,
	       stackBase<State> * _ldclist = NULL,
	       int nsteps = 0,
	       int nsnaps = 0)
      : step(_step), term(_term), init(_init), ldclist(_ldclist), _steps(nsteps),
       _snaps(nsnaps){}
  

    ~StdSimData() {}
    TimeStep<State> & getStep() const { return step; }
    TimeTerm & getTerm() const { return term; }
    Algorithm & getInit() const { return init; }

    stackBase<State> * getStack() const { return ldclist;}
	
    int getNsteps() const { 
      if (_steps > 0) return _steps;
      else {
	RVLException e;
	e<<"\n called from StdSimData::getNsteps(), nsteps ="<<_steps<<" (not positive) \n";
	throw e;
      }
    }

    int getNsnaps() const {return _snaps;}
  };

  /** Basic forgetful simulator - essentially a LoopAlg. */
  template<typename State>
  class StdSim: public Sim<State> {

  private:

    StdSim();

  public: 

    /** main constructor */
    StdSim(TimeStep<State> & step,
	   TimeTerm & term, 
	   Algorithm &initstep)
      : Sim<State>(step, term, initstep) {}
    
    /** use arg struct - for generic policy */
    StdSim(StdSimData<State> const & d)
      : Sim<State>(d.getStep(),d.getTerm(),d.getInit()) {}

    StdSim(StdSim<State> const & s)
      : Sim<State>(s) {}
    
    virtual ~StdSim() {}
    
    void run() {
      try {	

	LoopAlg a(this->step, this->term);
	ListAlg aa(this->initstep, a);
	aa.run();    
      }
      catch (RVLException & e) {
	e<<"\ncalled from StdSim::run\n";
	throw e;
      }

    }

     
  };



  /** Modified forgetful simulator, for adaptive time stepping. 
      Uses the class peekTimeStep to look into final time, 
      which could be used to alter the timestep length 
      to stop at the time dictated by terminator object
  */
 
  template<typename State>
  class StdIntSim: public Sim<State> {

  private:

    StdIntSim();

  public: 

    // main constructor 
    StdIntSim(PeekTimeStep<State> & step,
	      TimeTerm & term, 
	      Algorithm &initstep)
      : Sim<State>(step, term, initstep) {}
    
    StdIntSim(StdIntSim<State> const & d)
      : Sim<State>(d.getStep(),d.getTerm(),d.getInit()) {}

    
    virtual ~StdIntSim() {}
    
    void run() {
      try {	
	
	bool result = true;
	this->initstep.run();
	while ( (!this->term.query()) && result ) {
	  
	  // pass target time to make sure no overstepping takes place
	  // (timestep class is responsible for using the targettime to stop
	  //  over-simulating)
	  dynamic_cast<PeekTimeStep<State>  &>(this->step).peekTT(this->term.getTargetTime());
	  result = this->step.run();
	  
	}
      }

      catch (RVLException & e) {
	e<<"\ncalled from StdSim::run\n";
	throw e;
      }

    }
    
  };

 

  /** Remember-All Simulator. myStack handles how states are saved */
  template<typename State>
  class RASim: public Sim<State> {

  private:

    RASim();

    /** stateHist is a pre-constructed subclass of stackBase. This avoids the
	issue of accommodating various arguments that a stackBase subclass might
	need in order to construct a real object */
    stackBase<State> & stateHist;

    bool filledSH;
 
    int nt; // # of statehist elements (minus the initial state)



  public: 

    RASim(TimeStep<State> & step,
	  TimeTerm & term,
	  Algorithm & initStep,
	  stackBase<State> & stateHist_)
      : Sim<State>(step,term, initStep), filledSH(false),  nt(0), stateHist(stateHist_)
    { }
    
 
    RASim(RASim<State> const & s)
      : Sim<State>(s), filledSH(s.filledSH), nt(s.nt), stateHist(s.stateHist)
    {}
    
    virtual ~RASim() {
      /*
      cerr<<"in RASim destructor\n";
      while (stateHist.size() > 0) {
	stateHist.pop_back();
      }
      */
    }


    /** macro that runs the simulation and saves states onto the stackbase object */
    void fillSH() {
      // run initialization step
      this->initstep.run();

      // push initial state to stack
      stateHist.push_back(this->step.getState());

      // reset counter
      nt = 0;

      #ifdef DEBUG

      this->step.getTime().write(cout);
      this->term.getTargetTime().write(cout);
      cout << "< : " << (this->step.getTime() < this->term.getTargetTime()) << endl;
      cout << "== : " << (this->step.getTime() == this->term.getTargetTime()) << endl;
      #endif
 

      //      cout << "step.time = " << endl;
      //      this->step.getTime().write(cout);
      //      cout << "Term.ttime = " << endl;
      //      this->term.getTargetTime().write(cout);

      
      while ( this->step.getTime() < this->term.getTargetTime() ) {
	
       #ifdef DEBUG

	cout << "< : " << (this->step.getTime() < this->term.getTargetTime()) << endl;
	cout << "== : " << (this->step.getTime() == this->term.getTargetTime()) << endl;

      #endif

	this->step.run();
	stateHist.push_back(this->step.getState());
	++nt;
      }

#ifdef DEBUG      
      cout << "stateHist.size() = " << stateHist.size() << endl;
#endif
      filledSH = true;

     }

    /** to reset the RASim, the stack has to be cleared first! 
	Function fillSH() takes care of the initstep.run() part */
    void resetSim() {
      stateHist.clear();
      filledSH = false;
      nt = 0;
    }


    void run() {
  
      if (!filledSH) {  fillSH(); }

      while (( this->step.getTime() > this->term.getTargetTime() )) {
	
	try {
	  this->step.setTime(this->term.getTargetTime());	
	  
	  State & r = this->step.getState();
	  r.copy(stateHist.at(nt));
	
	  // naturally, you want the next to last element when you call run
	  --nt; 
	  stateHist.pop_back(); 

	  //	  cout << "nt = " << nt << endl;
	  //	  cout << "sh.size() = " << stateHist.size() << endl;
	
	  if (nt == 0) {
	    //	    cout << "TRUE" << endl;
	    //	    this->step.getTime().write(cout);
	    resetSim();
	  }
    
	 
	}
	catch (RVLException & e) {
	  e<<"\ncalled from RASim::run\n";
	  throw e;
	}
      }
    }

  };
  
  
  /** Run-from-Current-State Simulator. Runs simulation from its
      current state, no reinitialization is performed */
  template<typename State>
  class RCSim: public Sim<State> {

  private:

    RCSim();

  public: 

    RCSim(TimeStep<State> & step,
	  TimeTerm & term, 
	  Algorithm &initstep)
      : Sim<State>(step, term, initstep) {}
    
    RCSim(RCSim<State> const & s)
      : Sim<State>(s) {}
    
    virtual ~RCSim() {}
    
    void run() {
      try {
	LoopAlg a(this->step, this->term);	
	a.run();       
      }
      catch (RVLException & e) {
	e<<"\ncalled from StdSim::run\n";
	throw e;
      }
    }
    
  };
  
  /** Basic Run-from-Current-State Simulator */
  template<typename State>
  class StdRCSim: public Sim<State> {

  private:

    StdRCSim();

  public: 

    /** main constructor */
    StdRCSim(TimeStep<State> & step,
	   TimeTerm & term, 
	   Algorithm &initstep)
      : Sim<State>(step, term, initstep) {}
    
    /** use arg struct - for generic policy */
    StdRCSim(StdSimData<State> const & d)
      : Sim<State>(d.getStep(),d.getTerm(),d.getInit()) {}

    StdRCSim(StdRCSim<State> const & s)
      : Sim<State>(s) {}
    
    virtual ~StdRCSim() {}

    void run() {
      try {      
	LoopAlg a(this->step, this->term);	
	a.run();       
      }
      catch (RVLException & e) {
	e<<"\ncalled from StdRCSim::run\n";
	throw e;
      }
    }
      
  };

  /** runs member Sim (reference simulator) until it reaches the
      current time of the member TimeStep. Since only Time info
      transferred, reference simulator and timestep objects may have
      different state types */
  template<typename State, typename RefState>
  class TargetCurrTimeStep: public Algorithm {
  private:
    Sim<RefState> & ref;
    TimeStep<State> const & ts;
    
    TargetCurrTimeStep();
    
  public:
    TargetCurrTimeStep(Sim<RefState> & _ref, TimeStep<State> const & _ts)
      : ref(_ref), ts(_ts) {
          
    }
    TargetCurrTimeStep(TargetCurrTimeStep<State,RefState> const & t)
      : ref(t.ref), ts(t.ts){}
    ~TargetCurrTimeStep() {}

    void run() {
      try {
	ref.setTargetTime(ts.getTime());
	ref.run();
      }
      catch (RVLException & e) {
	e<<"\ncalled from TargetCurrTimeStep::run\n";
	throw e;
      }
    }
  };

  /** runs member Sim until it reaches the next time of the member
      TimeStep. Since only Time info transferred, reference simulator
      and timestep objects may have different state types */
  template<typename State, typename RefState>
  class TargetNextTimeStep: public Algorithm {
  private:
    Sim<RefState> & ref;
    TimeStep<State> const & ts;
       
    TargetNextTimeStep();
    
  public:
    TargetNextTimeStep(Sim<RefState> & _ref, TimeStep<State> const & _ts)
      : ref(_ref), ts(_ts) {}
    TargetNextTimeStep(TargetNextTimeStep<State,RefState> const & t)
      : ref(t.ref), ts(t.ts) {}
    ~TargetNextTimeStep() {}

    void run() {
      try {

	//	cout << "ts.getTime()" << endl;
	//	ts.getTime().write(cout);

	//	cout << "ts.getNextTime()" << endl;
	//	ts.getNextTime().write(cout);

	//	ts.write(cout);

	ref.setTargetTime(ts.getNextTime());
	ref.run();		
      }
      catch (RVLException & e) {
	e<<"\ncalled from TargetNextTimeStep::run\n";
	throw e;
      }
    }
  };

}


#endif	
