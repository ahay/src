#ifndef __TEST_DYN1
#define __TEST_DYN1

#include "rn.hh"
#include "rns.hh"
#include "rt.hh"
#include "alg.hh"
#include "dtstep.hh"
#include "contentpackage.hh"
#include "local.hh"
     



namespace TSOpt {

  using RVL::RVLException;
  using RVL::ContentPackage;
  using RVL::LocalDataContainer;
  using RVL::UnaryLocalFunctionObject;
  using RVL::UnaryLocalFunctionObjectScalarRedn;

  /** \page tests TSOpt Test Suite
      
  Tests based on logistic equation \f[ \frac{du}{dt} = 1 - u^2, u(0) =
  0.  \f] discretized by Euler's method, with \f$\Delta t =
  0.1\f$. Produces u(k) \f$\simeq u(k\Delta t)\f$, k=0,...,10.
  
  Source files are testk.cc, for various k. Also hard.cc is hardcoded
  source for a command to generate a list of values u(k), i=0,...,10.
  
  The State class in all of these examples is RnState, a simple
  array-based class with the added internal data necessary to
  implement the time-oriented functions of a State.  

  The tests are set up to run regression-style, with results checked
  against stored reference results (files are in Orig, not removed by
  make clean). Most of the results (all of which are human-readable)
  will be stable against platform/compiler change, though some may
  change in the last digits due to change in roundoff and thus produce
  false negatives. The results of the regression is summarized in
  regress.rpt.

  As for other regression tests implemented in the TRIP build system, the command

  make regress

  builds the binaries and runs them to build the regression results.

  <ol> 
  <li> test1:
  tests assignment operator for StdDiscreteTime.</li> 

  <li> test1a: tests
  a constructor for the Dyn1 TimeStep class and correct functioning of
  the getTime() and getNextTime() methods.</li> 

  <li> test2: constructs
  Dyn1 and a TimeTerm (\ref Terminator) class for detecting the end of
  the simulation, builds a \ref LoopAlg, and runs the simulation.</li>
  
  <li> test3: same as 2, but constructs a StdSim instead of a \ref
  LoopAlg.</li> 

  <li> test4: same as 3, but constructs a StdSim1Jet -
  this requires TimeSteps for linearized and adjoint time steps,
  together with appropriate terminators.

  <li> test5: same as 4, but runs linearized
  sim instead.</li> 

  <li> test5a: same as 4, but runs the adjoint linearized sim
  instead. This is the first of four tests of the adjoint code: this
  first test uses the "recompute from initial data" random access
  policy, which is very expensive for all but small problems.</li>

  <li>test5a_ra: same as test5a, but uses the "remember-all" policy,
  that is, stores all time levels of the reference state in an array
  and accesses them as needed.</li>

  <li>test5a_cp: same as test5a, but uses the optimal checkpointing
  policy (Griewank, Appl. Math. Opt. v1, 1992). Makes little
  difference for small problems, but much less storage for modest
  additional computation in general.</li>

  <li>test5a_acp: same as test5a_cp, but uses an adaptive variant of
  optimal checkpointing suitable for adaptive time stepping (which
  this expl does not have!). Method explained in Marco Enriquez'
  thesis.</li>

  <li> test6: All previous tests used
  the same Dyn1 constructor, which takes a \ref
  \ref RVL::ContentPackage<float,int> (i.e. a simple array) as input, and
  similarly for the linearized TimeStep class D_Dyn1. This test uses
  instead functions which transfer data between the RnState of the
  simulator, and an \ref RVL::Vector (member of \ref RVL::RnSpace)
  object. These functions are implemented as policy classes.</li>
  
  <li> test7: Performs the adjoint test on the simulation operator and
  its adjoint. Calculates the difference between <Ax,y> -
  <x,A*y>. Uses bonehead "recompute from initial data" policy for
  adjoint state.  </li>

  <li> test8: Derivative test. Checks the quality of the derivative 
  (linearized) simulation by comparing it to the FD derivative </li>

  </ol>
  */
	
  /** TimeStep class for logistic equation, Euler's method */
  class Dyn1: public StdDFTStep<RnState> {

  private:

    ostream & fpout;
    float dt;
    bool verbose;

    void dyn1fwd(rn const & u, float dt);
    
    Dyn1();

  public:

    Dyn1(ContentPackage<float,size_t> const & _c, 
	 float _dt, 
	 bool _verbose=false, 
	 ostream & _fpout=cout);

    Dyn1(float _dt, 
	 bool _verbose=false, 
	 ostream & _fpout=cout);
    Dyn1(Dyn1 const & t);

    ~Dyn1() {}
    
    void run();

    ostream & write(ostream & str) const;

  
    
  };

  /** TimeStep class for linearized logistic equation */

  class D_Dyn1: public StdDFTStep<RnState> {

  private:

    ostream & fpout;
    float dt;
    bool verbose;

    // reference state - mutable
    RnState & ref;

    void dyn1lin(rn const & du, rn const & u, float dt);

  public:

    D_Dyn1(ContentPackage<float,size_t> const & _dc, 
	   RnState & _ref, 
	   float _dt, 
	   bool _verbose=false, 
	   ostream & _fpout=cout);

    D_Dyn1(RnState & _ref, 
	   float _dt, 
	   bool _verbose=false, 
	   ostream & _fpout=cout);
    D_Dyn1(D_Dyn1 const & t);

    ~D_Dyn1() {}
    
    void run();
    ostream & write(ostream & str) const;
  };


  /** TimeStep class for adjoint logistic equation */

  class A_Dyn1: public StdDBTStep<RnState> {

  private:

    ostream & fpout;
    float dt;
    bool verbose;

    // reference state - mutable
    RnState & ref;

    void dyn1adj(rn const & au, rn const & u, float dt);

  public:

    A_Dyn1(ContentPackage<float,size_t> const & _ac, 
	   RnState & _ref, 
	   float _dt, 
	   bool _verbose=false, 
	   ostream & _fpout=cout);

    A_Dyn1(RnState & _ref, 
	   float _dt, 
	   bool _verbose=false, 
	   ostream & _fpout=cout);

    A_Dyn1(A_Dyn1 const & t);

    ~A_Dyn1() {}
    
    void run();
    ostream & write(ostream & str) const;
  };

  class NullAlg : public Algorithm {  
   
  public: 
    void run() {} 
  };


  template <typename T>
  class NullStateAlg : public StateAlg<T> {
  private:
    T & state;
     
  public:
    NullStateAlg(T & init): state(init) {}
    ~NullStateAlg() {}
    
    // void setState(const T & x) { state = x;}
    const T & getState() const {return state;}
    T & getState() {return state;}
    void run() {}

  };   

  template <typename T>
  class initFwd_GSL : public StateAlg<T> 
  {
  private:
    T _initState;
    StdTimeStep<T> & _d;
    
  public:
    initFwd_GSL(StdTimeStep<T> & d) :   _d(d)
    { 
      
      if (_d.getState().getrealrn().u) {_initState.copy(_d.getState()) ; }
      else {_initState.initialize(1,1,0.0);}
    }

    initFwd_GSL (initFwd_GSL const & s) : _initState(s._initState), _d(s._d) { }
    ~initFwd_GSL(){}

    // void setState(const T & x) { _initState = x;}
    const T & getState() const {return _initState;}
    T & getState() {return _initState;}

    void run() {
      try {
	RealTime t0;
	t0 = 0.0;

	//	cout << "~~~~~~ INIT STATE" << endl;
	//        _initState.write(cout);

	T & dummy(_d.getState());
	dummy.copy( _initState );

	//	cout << "init state after copy " << endl;
	//	_d.getState().write(cout);
	
	_d.setTime(t0);
      }
      catch (RVLException & e) {
	e<<"\ncalled from initFwd::run() \n";
	throw e;
      }
    }
    
  };

	
 template <typename T>
  class initFwd : public StateAlg<T> 
  {
  private:
    T _initState;
    Dyn1 & _d;
        
  public:
    initFwd(Dyn1 & d) :   _d(d), _initState(d.getState())
    { 
      //      _initState.copy(_d.getState()) ; 
      
    }
    initFwd (initFwd const & s) : _initState(s._initState), _d(s._d) { }
    ~initFwd(){}

    // void setState(const T & x) { _initState = x;}
    const T & getState() const {return _initState;}
    T & getState() {return _initState;}

    void run() {
      try {
	StdDiscreteTime t0;
	t0=0;
	//_d.setState(_initState);

	//cout << "~~~~~~ INIT STATE" << endl;
        //_initState.write(cout);

	T & dummy(_d.getState());
	dummy.copy( _initState );

	//_d.getState() = _initState;
	// _d.getState().setTime(t0);
      }
      catch (RVLException & e) {
	e<<"\ncalled from initFwd::run() \n";
	throw e;
      }
    }
    
  };


  template <typename T>
  class initLin : public StateAlg<T> 
  {
  private:
    T _initState;
    D_Dyn1 & _d;


  public:
    initLin(D_Dyn1 & d) :  _d(d), _initState(d.getState()) {}
			   //    { _initState.copy(d.getState());}
    initLin (initLin const & s) : _initState(s._initState), _d(s._d) { }
    ~initLin(){}

    //void setState(const T & x) { _initState = x;}
    const T & getState() const {return _initState;}
    T & getState() {return _initState;}

    void run() {
      try {
		cout << "@@@@@@ initLin" << endl;
	StdDiscreteTime t0;
	t0=0;
	//_d.setState(_initState);
	
	T & dummy(_d.getState());
	dummy.copy( _initState );

	//_d.getState()=_initState;
	//_d.getState().setTime(t0);

      }
      catch (RVLException & e) {
	e<<"\ncalled from initLin::run() \n";
	throw e;
      }
    }
    
  };



  template <typename T>
  class initAdj : public StateAlg<T> 
  {
  private:
    T _initState;
    A_Dyn1 & _d;
    

  public:
    initAdj(A_Dyn1 & d) :  _d(d), _initState(d.getState()) { 
      /*
      cerr<<"*******************************************"<<cerr;
      cerr<<"initAdj::_initState\n"; 
      _initState.write(cerr);
      cerr<<"initAdj::d\n"; 
      d.write(cerr);
      cerr<<"*******************************************"<<cerr;
      */
    }
			   //    { _initState.copy(d.getState());}
    initAdj (initAdj const & s) : _initState(s._initState), _d(s._d) { }
    ~initAdj(){}

    // void setState(const T & x) { _initState = x;}
    const T & getState() const {return _initState;}
    T & getState() {return _initState;}

    void run() {
      try {

	StdDiscreteTime t0;
	t0=9;
	//_d.setState(_initState);
	_d.getState().copy(_initState);
	_d.getState().setTime(t0);

      }
      catch (RVLException & e) {
	e<<"\ncalled from initAdj::run() \n";
	throw e;
      }
    }
    
  };

  template <typename T>
  class initFwd_RVL : public StateAlg<T>{
  private:
    T _initState;
    Dyn1 & _d;
   
    //    T & getState() {}
    //    const T & getState() const {}

  public:
    initFwd_RVL(Dyn1 & d) :  _d(d)
    {
      _initState.initialize(1,1,0);
      
    }
    initFwd_RVL (initFwd_RVL const & s) : _initState(s._initState), _d(s._d) { }
    ~initFwd_RVL(){}

    const T & getState() const {return _initState;}
    T & getState() {return _initState;}

    void run() {
      try {

	// cout << "CALLING INITFWD_RVL " << endl;
	
	StdDiscreteTime t0;
	t0=0;

	T & dummy(_d.getState());
	dummy.copy( _initState );

	
      }
      catch (RVLException & e) {
	e<<"\ncalled from initFwd_RVL::run() \n";
	throw e;
      }
    }
    
  };

  template <typename T>
  class initLin_RVL : public StateAlg<T>{
  private:
    T _initState;
    D_Dyn1 & _d;
    

  public:
    initLin_RVL(D_Dyn1 & d) :  _d(d)
    { _initState.initialize(1,1,0);}
    initLin_RVL (initLin_RVL const & s) : _initState(s._initState), _d(s._d) { }
    ~initLin_RVL(){}

    const T & getState() const {return _initState;}
    T & getState() {return _initState;}

    void run() {
      try {

	StdDiscreteTime t0;
	t0=0;
	_d.getState() = _initState;
	_d.getState().setTime(t0);
	//	_d.setState(_initState);

      }
      catch (RVLException & e) {
	e<<"\ncalled from initLin_RVL::run() \n";
	throw e;
      }
    }
    
  };


  template <typename T>
  class initAdj_RVL : public StateAlg<T>{
  private:
    T _initState;
    A_Dyn1 & _d;
    

  public:
    initAdj_RVL(A_Dyn1 & d) :  _d(d)
    { _initState.initialize(1,1,0);}
    initAdj_RVL (initAdj_RVL const & s) : _initState(s._initState), _d(s._d) { }
    ~initAdj_RVL(){}

    const T & getState() const {return _initState;}
    T & getState() {return _initState;}

    void run() {
      try {

	StdDiscreteTime t0;
	t0=9;
	_d.getState() = _initState;
	_d.getState().setTime(t0);
	// _d.setState(_initState);

      }
      catch (RVLException & e) {
	e<<"\ncalled from initAdj_RVL::run() \n";
	throw e;
      }
     
    }
    
  };







}
  
#endif
