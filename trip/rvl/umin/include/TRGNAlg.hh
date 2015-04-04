#ifndef RVL_UMIN_TRALG_
#define RVL_UMIN_TRALG_

#include "table.hh"
#include "op.hh"
#include "alg.hh"
#include "ioterm.hh"

namespace RVLUmin {

  using namespace RVL;
  using namespace RVLAlg;

  /** Generic trust region (truncated GN) step.

      Trust region decision parameters are constants, so pass by value.
      
      Operator evaluation, Function value, absolute and relative
      gradient norms, predicted and actual reduction, and trust radius
      are all altered by this object, so it stores mutable references
      for these items.

      Parameters / data members:
      <ul>
      <li> pol - solver policy, see description in TRGNAlg docs (TRGNAlg is a policy subclass). Generates LS solver via build method. On construction, an LS solver is required to initialize the residual norm and normal residual norm, passed by reference, as well as the RVL::Vector update.</li>
      <li> opeval - RVL::OperatorEvaluation, mutable reference, 0.5*getValue().normsq() = objective function</li>
      <li> eta1 - lower G-A parameter; when actred/predred < eta1, trust radius reduced. Typical value = 0.01</li>
      <li> eta2 - upper G-A parameter; when actred/predred > eta2 and trust region constraint is binding, trust radius increased. Typical value = 0.9</li>
      <li> gamma1 - trust region reduction factor. Typical value = 0.5</li>
      <li> gamma2 - trust region increase factor. Typical value = 1.8</li>
      <li> predred - predicted reduction - mutable reference</li>
      <li> actred - actual reduction - mutable reference</li>
      <li> jval - objective function value - mutable reference</li>
      <li> agnrm - gradient norm - mutable reference</li>
      <li> rgnrm - gradient norm scaled by reciprocal of initial (on constr) - mutable reference</li>
      <li> Delta - trust radius - mutable reference</li>
      <li> str - verbose output stream</li>
      </ul>
      
      Outline of run() method:
      <ul>
      <li> call trust region least squares solver to return step</li>
      <li> create trial iterate (storing input iterate in case of rejection)</li>
      <li> compute predicted, actual reductions in objective (predred, actred)</li>
      <li> if actred/predred < eta1, reduce trust radius by factor
      gamma1, reject step, restore iterate from stored input </li>
      <li> if actred/predred >= eta1, accept step (discard stored input)
      <ul>
      <li> if additionally actred/predred > eta2 and trust region constraint was binding on LS solution, increase trust radius by factor gamma2</li>
      </ul>
      <li> return</li>
      </ul>

  */
  template<typename Scalar, typename Policy >
  class TRGNStep: public Algorithm {
    
    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;
    
  public:
    
    /** here op = b-F(x), i.e. the residual - use ResidualOp wrapper 

	NOTE: this is the "optimistic" no-copy form - i.e. assumes
	that the step is likely to succeed, so the
	restor-xsave-compute branch is unlikely, hence not costly. A
	"pessimistic" form would be a mirror-image. The optimal form
	would involve a deep copy operation on OpEvals, which is
	certainly possible future mod.

	the quadratic model is
	
	0.5*\|b-F(x)\|^2 + p^Tg(x) + 0.5p^TH(x)p

	where H(x)=DF(x)^TDF(x), and H(x)p = g(x) = DF(x)^T(b-F(x)) (approx)

	The predicted reduction is

	0.5*\|b-F(x)\|^2 - 0.5p^Tg(x)
    */
    void run() {
      //      cerr<<"TRGNStep::run\n";
      // run then delete solver - first check whether trust region was
      // transgressed
      solver->run();  
      bool active = false;
      {
	Terminator & tsolver = dynamic_cast<Terminator &>(*solver);
	active = tsolver.query();
      }
      delete solver;	
      solver=NULL;

      // test for nontrivial step - failure is fatal
      if (p.norm()<minstep*pol.Delta) {
	nostep=true;
	return;
      }

      // here rnorm is final residual from TRLS solver
      atype jvalsave = jval;
      predred = jvalsave - 0.5*rnorm*rnorm;
      str<<"** predicted reduction = "<<predred<<" active = "<<active<<"\n";
      str<<"** length of step = "<<p.norm()<<"\n";
      //	cerr<<"TRGNStep: predred="<<predred<<endl;
      // provisionally update eval point - save current point in case 
      // need to revert
      Vector<Scalar> & x = opeval.getPoint();
      xsave.copy(x);
      x.linComb(-ScalarFieldTraits<Scalar>::One(),p);
      // reset (reconstruct) solver
      //      cerr<<"TRGNStep::run - main callg before call to pol.build Delta="<<pol.Delta<<endl;
      //      solver = pol.build(p,opeval.getDeriv(),opeval.getValue(),rnorm,nrnorm,str);	
      // by convention this should have recalculated the objective
      // function, gradient norm
      //      jval = 0.5*rnorm*rnorm;
      jval = 0.5*opeval.getValue().normsq();

      // compute actual reduction
      actred = jvalsave - jval;
      // Case 1: fail to satisfy G-A : rescind eval point update, shrink Delta
      str<<"** actual reduction: old jval = "<<jvalsave<<" new jval = "<<jval<<" actred = "<<actred<<"\n";
      if (actred < eta1*predred) {
	pol.Delta *= gamma1;
	x.copy(xsave);
      }
      // Case 2: otherwise keep update
      // Case 2.1: only grow Delta if step hit TR boundary
      if ((actred > eta2*predred) && active) { 
	int dcount = 0;
	Vector<Scalar> plast(opeval.getDomain());
	while (active && dcount<5) {
	  str<<"active trust region constraint and good step - attempt to increase TR\n";
	  atype deltalast=pol.Delta;
	  atype jvallast = jval;
	  
	  // save last update vector
	  plast.copy(p);
	  // stretch it
	  p.scale(gamma2);
	  // stretch TR
	  pol.Delta *= gamma2; 
	  str<<"last TR = "<<deltalast<<" new TR = "<<pol.Delta<<"\n";
	  // restore x
	  x.copy(xsave);
	  // compute DF(x)p - x not updated
	  Vector<Scalar> dp(opeval.getRange());
	  opeval.getDeriv().applyOp(p,dp);
	  // create quadratic residual
	  dp.linComb(-ScalarFieldTraits<Scalar>::One(),opeval.getValue());
	  // estimate reduction - note jvalsave is resid at x
	  predred = jvalsave - 0.5*dp.normsq();
	  str<<"new predred="<<predred<<"\n";
	  // update x
	  x.linComb(-ScalarFieldTraits<Scalar>::One(),p);
	  // create real residual
	  jval = 0.5*opeval.getValue().normsq();
	  str<<"last jval = "<<jvallast<<" new jval = "<<jval<<"\n";
	  // compute actual reduction
	  actred = jvalsave - jval;
	  str<<"new actred = "<<actred<<"\n";
	  // if G-A fails, back off
	  if (actred < eta1*predred) {
	    str<<"G-A failed, revert to previous estimate, unset active flag\n";
	    pol.Delta = deltalast;
	    jval=jvallast;
	    p.copy(plast);
	    x.copy(xsave);
	    x.linComb(-ScalarFieldTraits<Scalar>::One(),p);
	    active=false;
	  }	  
	  else {
	    str<<"G-A passed, try it again\n";
	  }
	  dcount++;
	}
      }

      // rejigger everything for the new point
      solver = pol.build(p,opeval.getDeriv(),opeval.getValue(),rnorm,nrnorm,str);		  
      agnrm=nrnorm;
      rgnrm=nrnorm*gnrmrecip;

      //      cerr<<"at end of TRGNStep::run - Delta = "<<pol.Delta<<endl;
    }
    
    atype getDelta() { return pol.Delta; }

    /** constructor 

	Depends on least-squares solver, chosen by the policy
	argument. Type of this object passed as template
	parameter. For characteristics, see TRGNAlg docs.

	Pass already-initialized operator evaluation to this
	method. This evaluation stores the RHS (opeval.getValue()) and
	LinearOp (opeval.getDeriv) of the G-N least squares
	problem. The evaluation point is updated as part of the step,
	so the opeval reference must be mutable.

	@param _pol - solver policy, see descriptions in TRGNAlg docs
	@param _opeval - updated reference, 0.5*getValue().normsq() = objective function
	@param _eta1 lower G-A parameter > 0
	@param _eta2 upper G-A parameter > eta1
	@param _gamma1 trust region reduction factor < 1
	@param _gamma2 trust region expansion factor > 1, gamma1*gamma2 < 1 
	@param _minstep min permitted step as a fraction of trust radius
	@param _nostep true if step norm was too small rel trust radius
	@param _predred predicted reduction - updated reference
	@param _actred actual reduction - updated reference
	@param _jval objective function value - updated reference
	@param _agnrm gradient norm - updated reference
	@param _rgnrm gradient norm scaled by reciprocal of original (on constr) - updated reference
	@param _str verbose output stream

	Requirements on parameters:
	<ul>
	<li>Delta >= 0</li>
	<li>0 < gamma1 < 1 < gamma2 </li>
	<li>gamma1 * gamma2 < 1</li>
	<li>0 < eta1 < eta2</li>
	<li>opeval = operator defining nonlinear least squares problem, evaluated at current iterate. used to supply RHS (getValue) and LinearOp (getDeriv) of G-N problem on call, and to evaluate residual and normal residual at G-N solution on return.</li>
	<li>pol - TRLS solver policy, as described above</li>
	</ul>

    */
    TRGNStep(Policy const & _pol,
	     OperatorEvaluation<Scalar> & _opeval,
	     atype _eta1,
	     atype _eta2,
	     atype _gamma1,
	     atype _gamma2,
	     atype _minstep,
	     bool & _nostep,
	     atype & _predred,
	     atype & _actred,
	     atype & _jval,
	     atype & _agnrm,
	     atype & _rgnrm,
	     ostream & _str) 
      : pol(_pol),
	opeval(_opeval), 
	eta1(_eta1),eta2(_eta2),gamma1(_gamma1),gamma2(_gamma2),minstep(_minstep),
	nostep(_nostep),predred(_predred), actred(_actred),
	jval(_jval), agnrm(_agnrm), rgnrm(_rgnrm),
	p(_opeval.getDomain()), xsave(_opeval.getDomain()), solver(NULL), str(_str) {
      
      // sanity test params
      if ( (gamma1 <= ScalarFieldTraits<atype>::Zero()) ||
	   (gamma1 > ScalarFieldTraits<atype>::One())  ||
	   (gamma2 < ScalarFieldTraits<atype>::One())  ||
	   (gamma1*gamma2 > ScalarFieldTraits<atype>::One()-100.0*numeric_limits<atype>::epsilon())  ||
	   (eta1  < ScalarFieldTraits<atype>::Zero())  ||
	   (eta1 > eta2) ) {
	RVLException e;
	e<<"Error: TRTRGNStep constructor\n";
	e<<"insane TR param inputs\n";
	e<<"gamma1= "<<gamma1<<"\n";
	e<<"gamma2= "<<gamma2<<"\n";
	e<<"eta1  = "<<eta1<<"\n";
	e<<"eta2  = "<<eta2<<"\n";
	throw e;
      }

      // no step taken yet so set reductions to zero
      actred = ScalarFieldTraits<atype>::Zero();
      predred = ScalarFieldTraits<atype>::Zero();
      
      // initial solver assignment
      solver = pol.build(p,opeval.getDeriv(),opeval.getValue(),rnorm,nrnorm,str);

      // solver policy: least squares solution of Ax=b. 
      // pass b=opeval.getValue, A=opeval.getDeriv references - cause reinitialization if needed
      // on creation, initialize x=0, r=b, rnorm=b.norm, nrnorm=A^Tb.norm 
      jval = 0.5*rnorm*rnorm;
      agnrm = nrnorm;

      // record initial reciprocal grad norm for future ref
      if (ProtectedDivision<atype>(ScalarFieldTraits<atype>::One(),nrnorm,gnrmrecip)) {
	// reset grad norms to zero
	agnrm = ScalarFieldTraits<atype>::Zero();
	rgnrm = ScalarFieldTraits<atype>::Zero();
	gnrmrecip = ScalarFieldTraits<atype>::Zero();
	str<<"NOTE: TRGNStep constructor\n";
	str<<"  initial value of normal residual norm = "<<nrnorm<<" below division threshold\n";
	str<<"  set residual norm values to zero, return - this should stop iteration\n"; 
      }
      else {
	rgnrm = ScalarFieldTraits<atype>::One();
      }
      // don't touch initial Delta

    }

    ~TRGNStep() { if (solver) delete solver; }

  private:
    
    Policy const & pol;
    OperatorEvaluation<Scalar> & opeval;
    atype & predred;
    atype & actred;
    atype & jval;
    atype & agnrm;
    atype & rgnrm;
    bool & nostep;
    atype gnrmrecip;
    atype eta1;
    atype eta2;
    atype gamma1;
    atype gamma2;
    atype minstep;
    atype rnorm;
    atype nrnorm;
    Vector<Scalar> p;
    Vector<Scalar> xsave;
    // workspace for solver
    mutable Algorithm * solver;
    ostream & str;
  };

  /** Trust Region iteration. Approximates local solution of \f[\min_x
      \|f(x)\|^2\f] by a variant of the Steihaug-Toint trust region
      globalization of the Gauss-Newton algorithm, with flexible
      specification of the least-squares substep. See for instance
      T. Steihaug, "The conjugate gradient method and trust regions in
      large scale optimization", SIAM Journal on Numerical Analysis,
      v. 20 (1983), pp. 626-637; A. Conn, N. Gould, and P. Toint,
      <i>Trust Region Methods</i>, SIAM, 2000; and J. Nocedal and
      S. Wright, <i>Numerical Optimization</i>, Spring, 1999, esp. Ch
      4.
      
      Note that this algorithm will also approximate solutions of
      \f[\min_x \|g(x)-b\|^2,\f] provided that \f$f(x)=g(x)-b\f$ is
      passed to the constructor.

      The iteration updates the current solution estimate \f$x_c\f$ to
      a new estimate \f$x_+\f$ by solving approximately the
      <i>trust-region Gauss-Newton subproblem</i> for the increment
      \f$\delta x\f$:

      \f[ \delta x = \mbox{argmin }\|Df(x_c)\delta x + f(x_c)\|^2
      \mbox{ subject to } \|\delta x \| \le \Delta \f]

      after which the update is computed as \f$x_+=x_c + \delta x\f$. 

      The <i>trust radius</i> represents a size estimate of the region
      near \f$x_c\f$ in which the Gauss-Newton quadratic is a good
      approximation to the least-squares objective. The genius of the
      algorithm lies in the rules for updating \f$\Delta\f$: it is
      decreased by a factor < 1 when the actual decrease in the
      objective obtained by \f$ x_c \mapsto x_+\f$ is substantially
      less than the decrease predicted by the Gauss-Newton model,
      increased by a factor > 1 when the actual decrease is close to
      the predicted decrease <i>and</i> the trust region constraint is
      binding, and otherwise left alone. The step algorithm
      RVLUmin::TRGNStep manages \f$\Delta\f$ - see its documentation
      for detailed information.
      
      TRGNAlg depends on a choice of least squares solver with trust
      region truncation. This choice is passed <i>by policy</i>,
      that is, both by inheritance (mixin) and as a template
      parameter. The policy type must define a public member function
      with following signature: <p> [subtype of Algorithm and
      Terminator] * build(...) const; <p> with arguments
      
      <ul> <li> Vector<Scalar> &, // on return, output of TR LS
      solve</li> <li> OperatorEvaluation<Scalar> &, // eval object -
      offers RHS and LinearOp of GN problem</li> <li> AbsScalar
      &, // on return, residual norm</li> <li> AbsScalar &, // on
      return, normal residual norm</li> <li> AbsScalar, // trust
      radius</li> <li> ostream & // verbose output unit</li> </ul>

      Policies should be default-instantiated, with no arguments,
      and should exhibit some sensible default behaviour (likely
      useless). For useful behaviour, some policies may require
      additional information - should be supplied in driver code by
      calling appropriate additional methods of particular policy
      classes.
      
      The least squares solver built by the policy must be supplied
      with the attributes of both RVLAlg::Algorithm and
      RVLAlg::Terminator. The Algorithm::run() method changes the
      states of the arguments as described above. The
      Terminator::query() method returns true if trust region
      truncation was applied, else false (i.e. least squares solution
      approximation was interior to the trust region).

      For example, the RVLUmin::CGNEPolicy class implements a policy
      as just described. Use it to specialize TRGNAlg to create a
      <i>Gauss-Newton-Krylov</i> or Steihaug-Toint algorithm, as
      described in the references mentioned earlier. The conjugate
      gradient iteration parameters must be supplied after
      instantiation for nontrivial execution, so the
      RVLUmin::CGNEPolicy includes a suitable method (assign(...))
      whichh must be called on the TRGNAlg object after
      construction. See the <a href="../../testsrc/testtrgn.cc">TRGN
      functional test source</a> for an explicit example of this use
      mode.

      The step algorithm RVLUmin::TRGNStep manages both the solution
      of the least squares problem and the trust radius (data member
      Delta) - see documentation on RVLUmin::TRGNStep for details.

      Parameters supplied to constructor. Absolute value type
      (eg. double for complex<double> denoted as abstype. See
      constructor docs for complete list:

      <ul>

      <li>x - RVL::Vector object, mutable reference. On construction,
      initial estimate of solution. On return, final estimate of
      solution returned by algorithm</li>

      <li>op - RVL::Operator object, const reference - operator
      defining the functional equation to be solved in the
      least-squares sense. As noted above, the functional equation
      takes the form \f$f(x)=0\f$, so any nontrivial "right-hand side"
      must be folded into the definition of \f$f\f$.</li>

      <li> _maxcount - int, maximum permitted TR-GN steps. Typical value = 10</li>

      <li> _jtol - abstype, stopping threshold for objective
      (i.e. one-half mean square of \f$f\f$). In some cases, it may be
      reasonable to stop when the objective gets small enough. To
      disable, set = zero. Typical value = 0.0</li>

      <li> _agtol - abstype, absolute stopping threshold for gradient
      norm, usable when scale information about the gradient is
      available, otherwise can be set to zero. Typical value = 0.0</li>

      <li> _rgtol - abstype, relative stopping threshold for gradient
      norm. Does not require scale of gradient or objective to be
      known, but may result in failure to converge if initial solution
      estimate is already accurate. Typical value = 0.01</li>

      <li> _initDelta - abstype, initial trust radius, requires the
      same kind of usually unavailable knowledge about the solution scale to
      set intelligently, as does the initial step in a line search
      method. Choose some convenient number and let the trust region
      algorithm take care of adjustments. Typical value = 1.0.</li>

      <li> _str - ostream, verbose output</li>

      <li> other parameters proper to trust region algorithm implemented
      in RVLUmin::TRGNStep - see its docs for details</li>

      </ul>

      Usage: instantiate; call post-construction initialization of
      least-squares solver policy parent class if necessary to assign
      additional attributes; call RVLUmin::TRGNAlg::run().

      Typical use case: see <a href="../../testsrc/testtrgn.cc">TRGN
      functional test source</a>
    
  */
  template<typename Scalar, typename Policy>
  class TRGNAlg: public Algorithm, public Policy {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  public:

    void run() {
      try {
	bool nostep = false;
	TRGNStep<Scalar, Policy> step(*this,opeval,eta1,eta2,gamma1,gamma2,minstep,
				      nostep,predred,actred,jval,agnrm,rgnrm,str);
	BoolTerminator bstop(nostep);
	VectorCountingThresholdIterationTable<atype> vstop(maxcount,names,nums,tols,str);
	vstop.init();
	atype jval0=jval;
	atype agnrm0=agnrm;
	atype delta0=step.getDelta();
	OrTerminator stop(bstop,vstop);
	LoopAlg doit(step,stop);
	doit.run();
	actcount=vstop.getCount();
	str<<"\n******************* TRGN Summary *******************\n";
	if (nostep && maxcount > 0) { 
	  str<<"TRGN: no update from GN step at TR step "<<actcount<<"\n";
	  str<<"  possibly gradient already below threshhold\n";
	}
	if (actcount > 0) {
	  str<<"iteration count            = "<<actcount<<"\n";
          str<<"initial objective          = "<<jval0<<"\n";
          str<<"final objective            = "<<jval<<"\n";
          str<<"objective redn             = "<<jval/jval0<<"\n";
          str<<"initial gradient norm      = "<<agnrm0<<"\n";
          str<<"gradient norm              = "<<agnrm<<"\n";
          str<<"gradient redn              = "<<rgnrm<<"\n";
	  str<<"initial trust radius       = "<<delta0<<"\n";
	  str<<"final trust radius         = "<<step.getDelta()<<"\n";
	}	  
	else {
	  if (!nostep) {
	    str<<"TRGN no initial step due to internal error\n";
	  }
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from TRGNAlg::run\n";
	throw e;
      }
    }
    
    /** convenience function for test codes */
    int getCount() const { return actcount; }
    
    /** constructor, first form.
	
	  Parameters:
	  @param op - operator figuring in least-squares problem
	  @param x - solution Vector - initial guess on call, estimated solution on return
	  @param _maxcount - max permitted TR steps (calls to TRGNStep)
	  @param _jtol - stopping tolerance for least squares residual
	  @param _agtol - absolute stopping tolerance for gradient norm
	  @param _rgtol - relative stopping tolerance for gradient norm
	  @param _eta1 - lower G-A parameter
	  @param _eta2 - upper G-A parameter
	  @param _gamma1 - trust region reduction factor
	  @param _gamma2 - trust region expansion factor
	  @param _initDelta - initial trust radius 
	  @param _str - verbose output stream
	
	  requirements on parameters:
	  <ul>
	  <li>Delta >= 0</li>
	  <li>0 < gamma1 < 1 < gamma2 </li>
	  <li>gamma1 * gamma2 < 1</li>
	  <li>0 < eta1 < eta2</li>
	  <li>opeval = operator defining nonlinear least squares problem, evaluated at current iterate. used to supply RHS (getValue) and LinearOp (getDeriv) of G-N problem on call, and to evaluate residual and normal residual at G-N solution on return.</li>
	  <li>pol - TRLS solver policy, as described above</li>
	  </ul>
    */
    TRGNAlg(Operator<Scalar> const & op, 
	    Vector<Scalar> & x,
	    int _maxcount,
	    atype _jtol,
	    atype _agtol,
	    atype _rgtol,
	    atype _eta1,
	    atype _eta2,
	    atype _gamma1,
	    atype _gamma2,
	    atype _minstep,
	    ostream & _str)
      : Policy(), opeval(op,x),
	jtol(_jtol),agtol(_agtol),rgtol(_rgtol),
	eta1(_eta1),eta2(_eta2),gamma1(_gamma1),gamma2(_gamma2),
        maxcount(_maxcount),minstep(_minstep),
	names(6),nums(6),tols(6),actcount(0),
	str(_str) { this->inittable(); }
    
    /** second form of constructor - passes most argument via a Table
	- see docs for first form of constructor for items which must
	be included in the Table and constraints upon them. */
    
    TRGNAlg(Operator<Scalar> const & op, 
	    Vector<Scalar> & x,
	    Table & t,
	    ostream & _str)
      : opeval(op,x),
	jtol(getValueFromTable<atype>(t,"ResidualTol")),
	agtol(getValueFromTable<atype>(t,"AbsGradTol")),
	rgtol(getValueFromTable<atype>(t,"RelGradTol")),
	eta1(getValueFromTable<Scalar>(t,"MinDecrease")),
	eta2(getValueFromTable<Scalar>(t,"GoodDecrease")),
	gamma1(getValueFromTable<Scalar>(t,"StepDecrFactor")),
	gamma2(getValueFromTable<Scalar>(t,"StepIncrFactor")),
	maxcount(getValueFromTable<int>(t,"MaxItn")), 
	minstep(getValueFromTable<int>(t,"MinStepTol")),
	names(6),nums(6),tols(6),actcount(0),
	str(_str) { 
      this->inittable(); }

  private:
    OperatorEvaluation<Scalar> opeval;
    mutable atype predred;
    mutable atype actred;
    mutable atype jval;
    atype jtol;
    mutable atype agnrm;
    mutable atype rgnrm;
    atype agtol;
    atype rgtol;
    atype eta1;
    atype eta2;
    atype gamma1;
    atype gamma2;
    int maxcount;
    atype minstep;
    vector<string> names;
    vector<atype *> nums;
    vector<atype> tols;
    mutable int actcount;
    ostream & str;

    void inittable() {
      names[0]="LS Objective"; nums[0]=&jval; tols[0]=jtol;
      names[1]="Gradient Norm"; nums[1]=&agnrm; tols[1]=agtol;
      names[2]="Rel Grad Norm"; nums[2]=&rgnrm; tols[2]=rgtol;
      names[3]="Actual Redn"; nums[3]=&actred; tols[3]=-numeric_limits<atype>::max();
      names[4]="Predicted Redn"; nums[4]=&predred; tols[4]=ScalarFieldTraits<atype>::Zero();
      names[5]="Trust Radius"; nums[5]=&(this->Delta); tols[5]=ScalarFieldTraits<atype>::Zero();
    }
  };
}
#endif
