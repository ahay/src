// chebalg.hh
// created by Yin 01/11/13
// last modified 01/25/13

// Much code is shamelessly stolen from the umin Cheb.H
// originally written by WWS, with his permission

/*************************************************************************

Copyright Rice University, 2004.
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


#ifndef __RVLALG_UMIN_Cheb_H
#define __RVLALG_UMIN_Cheb_H

#include "alg.hh"
#include "terminator.hh"
#include "linop.hh"
#include "table.hh"

using namespace RVLAlg;

namespace RVLUmin {

  using namespace RVL;
  using namespace RVLAlg;    

  template<typename Scalar>
  void ChebCoeff(typename ScalarFieldTraits<Scalar>::AbsType alpha, 
		 typename ScalarFieldTraits<Scalar>::AbsType epsilon, 
		 typename ScalarFieldTraits<Scalar>::AbsType gamma, 
		 int & kmax,
		 std::vector<typename ScalarFieldTraits<Scalar>::AbsType> & coeff) {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;
    // NOTE: estimating the necessary number of iterations
    atype one = ScalarFieldTraits<atype>::One();
    atype beta=0;
    atype epsest=0;
    ProtectedDivision<atype>(sqrt(alpha)*epsilon,one+sqrt(alpha),epsest);
    atype epsk;
    epsk= epsest + one;
    atype tmp = one;
    ProtectedDivision<atype>(one - gamma,one + gamma,beta);
    atype q=0;
    ProtectedDivision<atype>(beta,one + sqrt(one-beta*beta),q);
    int k = 1;
    while (epsk > epsest && k<=kmax) {
      tmp = tmp * q;
      ProtectedDivision<atype>(2*tmp,one + tmp*tmp,epsk);
      k = k + 1;
    }
    kmax = k-1;
    // NOTE: compute Chebyshev coefficients
    atype ckm = one;
    atype ck = one / beta;
    atype ckp=0;
    coeff.reserve(kmax+1);   // allocate memory for coeff
    coeff[0] = ScalarFieldTraits<atype>::Zero();
    coeff[1] = one;
    for (int i=2; i<kmax+1; i++) {
      ProtectedDivision<atype>(2*ck,beta,ckp);
      ckp = ckp - ckm;
      coeff[i] = one + ckm / ckp;
      ckm = ck;
      ck = ckp;
    }
  }

  /** Single step of Chebyshev iteration for a SPD op.
      
  On construction, internal workspace allocated and initialized.
  Each step updates internal state of ChebStep object. Since
  solution vector, residual norm, and normal residual norm are
  stored as mutable references to external objects, these external
  objects are updated as well.

  IMPORTANT NOTE: this version of the algorithm assumes that the
  solution vector reference (internal data member x) refers to a
  zero vector on initialization. To accommodate nontrivial initial
  guess, <i>modify the right-hand-side vector</i> (argument _b)
  externally.

  Solution vector (x), iteration count, residual norm, and normal
  residual (gradient) norm are all references to external objects,
  which may be monitored by appropriate terminators to build a
  LoopAlg out of this Algorithm.

  See ChebAlg for description of a fully functional algorithm
  class, combining this step with a Terminator to make a LoopAlg.
  */
  template<typename Scalar>
  class AltChebStep : public Algorithm {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  public:

    AltChebStep(LinearOp<Scalar> const & _A,
		Vector<Scalar> & _x,
		Vector<Scalar> const & _b,
		atype & _nrnorm,
		atype & _rhoest,    // input operator spectrum in case it is know prior
		std::vector<atype> &_coeff,
		int & _kc,
		atype _gamma = 0.04,     // inversion level
		atype _epsilon = 0.001,  // error reduction
		ostream & _str = cout)
      : A(_A), x(_x), b(_b), nrnorm(_nrnorm), coeff(_coeff), gamma(_gamma), epsilon(_epsilon), rhoest(_rhoest),kc(_kc), str(_str), dx(A.getDomain()), r(A.getDomain()), ndx(A.getDomain()) {

      kc=0;
      atype one = ScalarFieldTraits<atype>::One();
      atype tmp = one;

      if (b.getSpace() == A.getRange())
	A.applyAdjOp(b,r);
      else
	r.copy(b);
      nrnorm=r.norm();
      NormalLinearOp<Scalar> N(A);
      N.applyOp(r,ndx);

      tmp = abs(r.inner(ndx));
      if (ProtectedDivision<atype>(tmp,r.normsq(),RQ)) {
	RVLException e;
	e<<"Error: AltChebStep::AltChebStep() from ProtectedDivision: RQ\n";
	throw e;
      }
      // check spectral bound - if violated, reset. Note: this is NOT
      // a restart, as s has not been set yet!! Also will not be treated
      // as restart in alg, as it happens in the step constructor, not
      // in the run method!!

      // WWS 28.07.14: use Krylov-Bogoliubov upper bound rather than
      // ad-hoc scaling, here and in run method

      if(RQ > rhoest){
	Vector<Scalar> tmpvec(A.getDomain());
	tmpvec.copy(ndx);
	tmpvec.linComb(-RQ,r);
	atype tau;
	if (ProtectedDivision<atype>(tmpvec.norm(),nrnorm,tau)) {
	  RVLException e;
	  e<<"Error: AltChebStep::AltChebStep() from ProtectedDivision: K-B\n";
	  throw e;
	}
	rhoest = RQ + tau;
	str<<"Restart: estimated spectrum bound = "<< rhoest << endl;
      }


      // initial assignment of spectral scaling factor s
      ProtectedDivision<atype>(2*one,(one+gamma)*rhoest,s);

      x.zero();
      dx.zero();


      str<<"kc=0 nrnorm="<<nrnorm<<endl;
    }
      
    /**
       Run a single step of the AltChebyshev iteration for the normal equations
    */
    void run() {
      try {
	atype tmp;
	atype one = ScalarFieldTraits<atype>::One();

	dx.linComb(s*coeff[kc+1],r,coeff[kc+1]-one); // this = a*x+b*this
          
	x.linComb(one,dx);
	NormalLinearOp<Scalar> N(A);
	N.applyOp(dx,ndx);
          
	r.linComb(-1*one,ndx);
	nrnorm=r.norm();
	tmp = abs(dx.inner(ndx));
	if (ProtectedDivision<atype>(tmp,dx.normsq(),RQ)) {
	  RVLException e;
	  e<<"Error: AltChebStep::run() from ProtectedDivision: RQ\n";
	  throw e;
	} 
	kc++;
	str<<"kc="<<kc<<" nrnorm="<<nrnorm<<endl;
	if(RQ > rhoest){
	  Vector<Scalar> tmpvec(A.getDomain());
	  tmpvec.copy(ndx);
	  tmpvec.linComb(-RQ,dx);
	  atype tau;
	  if (ProtectedDivision<atype>(tmpvec.norm(),dx.norm(),tau)) {
	    RVLException e;
	    e<<"Error: AltChebStep::AltChebStep() from ProtectedDivision: K-B\n";
	    throw e;
	  }
	  rhoest = RQ + tau;
	  str<<"Restart: estimated spectral bound = "<< rhoest << endl;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from AltChebStep::run()\n";
	throw e;
      }
     
    }
    
    atype getSpectrumBound() const { return rhoest; }

    ~AltChebStep() {}

  private:

    // references to external objects
    LinearOp<Scalar> const & A;
    Vector<Scalar> & x;
    Vector<Scalar> const & b;
    atype & nrnorm;
    atype & rhoest;         // estimated spectrum bound
    int &kc;
    ostream & str;
    
    // added for AltChebyshev
    std::vector<atype> &coeff;
 
    atype gamma; 		   // inversion level
    atype epsilon; 		   // error reduction factor
    atype beta;
    atype RQ;
    atype s;

    // end of added
    
    // need four work vectors and one scalar as persistent object data
    Vector<Scalar> dx;
    Vector<Scalar> r;
    Vector<Scalar> ndx;
      
  };
  
 
  /** Chebyshev polynomial algorithm - efficient implementation for
      normal equations
      \f[ A^{\prime} A x = A^{\prime} b\f]
      for solving the linear least squares problem
      \f[ \min_{x} \vert A x - b \vert^2 \f].

      This is Chebyshev iteration Algorithm  as stated in Symes and
      Kern, Geophysical Prospecting vol. 42 pp. 565-614 1994 (see p. 578).
      We use variable names following Symes and Kern's notation insofar as
      possible. 

      Step 1:

      Choose an inversion level \f$\gamma\f$, typically 0.04;
      Choose an error reduction factor \f$\epsilon\f$;
      Choose a 'fudge factor' \f$\alpha>1.0\f$.

      Step 2: Compute the Chebyshev coefficients and estimate the 
      necessary number of iterations.

      \f$ \epsilon_{\mbox{est}}=\frac{\sqrt{\alpha}\epsilon}
      {1+\sqrt{\alpha}}\f$;

      \f$ \beta = \frac{1-\gamma}{1+\gamma}\f$;

      \f$ q = \frac{\beta}{1+\sqrt{1-\beta^2}} \f$.
        
      Define error reduction factor after \f$k\f$ step as
      \f$ \epsilon_k=\frac{2q^k}{1+q^{2k}}\f$.
      Let \f$ k_{\mbox{max}} \f$ be the smallest \f$k\f$ which
      satisfies \f$ \epsilon_k< \epsilon_{\mbox{est}} \f$.

      \f$ c_0 =1, c_1=\frac{1}{\beta},\omega_0=0,\omega_1=1. \f$
        
      For \f$k=1,\cdot,k_{\mbox{max}}-1\f$, compute
        
      \f$ c_{k+1} = \frac{2}{\beta}c_k-c_{k-1}\f$,
        
      \f$ \omega_{k+1}=1+\frac{c_{k-1}}{c_{k+1}}.\f$

      Step 3: Application of the Chebyshev polynomial by recursive
      application of the normal operator.
	
      a. Initialization
	
      \f$ x_0 = 0, dx_0 =0\f$,
	
      \f$ r_0 = A^{\prime}b, ndx_0=A^{\prime}Ae_0\f$.
        
      \f$ RQ_0 = \frac{\langle r_0, ndx_0\rangle }
      {\langle r_0, r_0 \rangle} \f$.
        
      \f$ \lambda_{\mbox{est}} = \alpha RQ_0, 
      s = \frac{2}{(1+\gamma)\lambda_{\mbox{est}}}\f$.
  
      b. Iteration
	
      For \f$k=0,\cdot, k_{\mbox{max}}-1\f$ 

      \f$ dx_{k+1} = (\omega_{k+1}-1)dx_k+s\omega_{k+1}r_k\f$,

      \f$ x_{k+1} = x_k + dx_{k+1}\f$,

      \f$ ndx_{k+1} = A^{\prime}A dx_{k+1}\f$,

      \f$ r_{k+1} = r_k-ndx_{k+1}\f$,

      \f$ RQ_{k+1}=\frac{\langle dx_{k+1},ndx_{k+1}\rangle}
      {\langle dx_{k+1},dx_{k+1}\rangle}\f$.

      if \f$ RQ_{k+1}> \lambda_{\mbox{est}}\f$, replace
      \f$\lambda_{\mbox{est}}\f$ by 
      \f$\alpha RQ_{k+1}\f$. Recompute 
      \f$ s = \frac{2}{(1+\gamma)\lambda_{\mbox{est}}}\f$,
      and restart step b.
	
      The final outputs are the estimated solution
      \f$x_{\mbox{est}}=x_{k_{\mbox{max}}}\f$
      and the estimated normal residual 
      \f$e_{\mbox{est}}=e_{k_{\mbox{max}}}\f$.

      Structure and function: Combines ChebStep with a Terminator
      which displays iteration count, residual norm, and normal
      residual norm on output stream (constructor argument _str), and
      terminates if iteration count exceeds max or residual norm or
      normal residual norm fall below threshhold (default =
      10*sqrt(macheps)). 

      Usage: construct ChebAlg object by supplying appropriate
      arguments to constructor. On return from constructor, solution
      vector initialized to zero, residual norm to norm of RHS, and
      normal residual norm to norm of image of RHS under adjoint of
      operator. Then call run() method. Progress of iteration written
      on output unit. On return from run(), solution vector stores
      final estimate of solution, and residual norm and normal
      residual norm scalars have corresponding values.

      Typical Use: see <a href="../../testsrc/testCheb.cc">
      functional test source</a>.

      IMPORTANT NOTE: This class is also an RVLAlg::Terminator
      subclass. 

      IMPORTANT NOTE: The solution vector and residual and normal
      residual scalars are external objects, for which this algorithm
      stores mutable references. These objects are updated by
      constructing a ChebAlg object, and by calling its run() method.

      IMPORTANT NOTE: this version of the algorithm initializes the
      solution vector to zero. To accommodate nontrivial initial
      guess, <i>modify the right-hand-side vector</i> (argument _rhs)
      externally.

      IMPORTANT NOTE: in order that this algorithm function properly
      for complex scalar types, a careful distinction is maintained
      between the main template parameter (Scalar) type and its
      absolute value type. All of the scalars appearing in the
      algorithm are actually of the latter type.

      See constructor documentation for description of parameters.


  */

  template<typename Scalar>
  class AltChebAlg: public Algorithm {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  public:

    /** Constructor
	
    @param _x - mutable reference to solution vector (external),
    initialized to zero vector on construction, estimated solution
    on return from ChebAlg::run().
    
    @param _inA - const reference to LinearOp (external) defining
    problem
    
    @param _rhs - const reference to RHS or target vector in
    either normal equations or basic system (i.e. can be either in
    domain or range of A) (external)
    
    @param _rnorm - mutable reference to residual norm scalar
    (external), initialized to norm of RHS on construction, norm
    of estimated residual at solution on return from
    ChebAlg::run()
    
    @param _nrnorm - mutable reference to normal residual (least
    squares gradient) norm scalar (external), initialized to morm
    of image of RHS under adjoint of problem LinearOp on
    construction, norm of estimated normal residual at solution on
    return from ChebAlg::run()
    
    @param _gamma - inversion level for Chebyshev
    
    @param _epsilon - stopping threshold for normal residual norm
    
    @param _alpha  - 'fudge factor'
    
    @param _maxcount - max number of iterations permitted, default
    value = 10
    
    default value = max absval scalar (which makes the trust
    region feature inactive)
    
    @param _str - output stream
    
    */
    AltChebAlg(RVL::Vector<Scalar> & _x, 
	       LinearOp<Scalar> const & _inA, 
	       Vector<Scalar> const & _rhs, 
	       atype & _nrnorm,
	       atype _gamma = 0.04,     // inversion level
	       atype _epsilon = 0.001,  // error reduction
	       atype _alpha = 1.1,      // 'fudge factor'
	       atype _rhoest=0.0,
	       int _maxcount = 10,      // upper bound of iterations
	       bool _restart = true,
	       ostream & _str = cout)  
      : inA(_inA), x(_x), rhs(_rhs), nrnorm(_nrnorm), gamma(_gamma),epsilon(_epsilon),alpha(_alpha),rhoest(_rhoest),maxcount(_maxcount), kmax(_maxcount), kc(0), ktot(0), restart(_restart), str(_str) {


      
      // NOTE: reference to coeff has been passed to ChebStep object step, however
      // coeff has not been initialized.
      ChebCoeff<Scalar>(alpha, epsilon, gamma, kmax, coeff); 
      cerr<<"kmax="<<kmax<<endl;
      for (int i=0;i<kmax;i++) {
	cerr<<"coeff["<<i<<"]="<<coeff[i]<<endl;
      }
    }
  

    ~AltChebAlg() {}
    
    void run() { 
      
      // need to put in alg form
      nrt=0;
      ktot=0;
      atype nrnorm0;
      kc=0;
      MaxTerminator<int> stop2(kc,kmax-1);
      while (!(stop2.query()) && ktot < maxcount) {
	// resets x=0, also possibly rhoest, computes nrnorm
	Algorithm * step = new AltChebStep<Scalar>(inA,x,rhs,nrnorm,rhoest,coeff,kc,gamma,epsilon,str);
	if (nrt==0) {
	  nrnorm0=nrnorm;
	}
	if (restart) {
	  atype curr_rhoest = rhoest;
	  MaxTerminator<atype> stop1(rhoest,curr_rhoest);
	  OrTerminator stop(stop1,stop2);
	  LoopAlg doit(*step,stop);
	  doit.run();
	}
	else {
	  LoopAlg doit(*step,stop2);
	  doit.run();
	}
	delete step;
	// success
	nrt++;
	ktot += kc;
      }
      if (!(stop2.query())) {
	// fatal error!
	RVLException e;
	e<<"Error: AltChebAlg::run\n";
	e<<"  failed to converge in allowed number of restarts\n";
	throw e;
      }
      if (ktot>=maxcount && !(stop2.query())) {
	// fatal error!
	RVLException e;
	e<<"Error: AltChebAlg::run\n";
	e<<"  failed to converge in allowed number of iterations\n";
	throw e;
      }    	  

      str<<"===== Restarted Chebyshev Iteration: Successful Completion =====\n";
      // display results
      str<<"\n ******* summary ********  "<<endl;
      str<<"initial gradient norm      = "<<nrnorm0<<endl;
      str<<"gradient norm              = "<<nrnorm<<endl;
      str<<"gradient redn              = "<<nrnorm/nrnorm0<<endl;
      str<<"final spec rad est         = "<<rhoest<<endl;
      str<<"restart flag               = "<<restart<<endl;
    }
    
    int getCount() const { return kc; }
    int getTotalCount() const { return ktot; }
    int getRestartCount() const { return nrt+1; }
    atype getSpectrumBound() const { return rhoest; }
    
    
  private:
    
    LinearOp<Scalar> const & inA;  // operator
    Vector<Scalar> & x;            // state - solution vector
    Vector<Scalar> const & rhs;    // reference to rhs
    atype & nrnorm;                // gradient norm
    
    // added for Chebyshev
    std::vector<atype> coeff;
    atype gamma; 		   // inversion level
    atype epsilon; 		   // error reduction factor
    atype alpha;		   // 'fudge factor'
    atype rhoest;                  // spectrum bound of op
    int maxcount;                  // upper bound for iteration count
    int nrt;                       // restarts
    bool restart;                  // flag
    
    /** Necessary number of iterations for a given error reduction factor.
	Special for Chebyshev iteration*/
    int kmax;
    int kc;
    int ktot;
    ostream & str;                 // stream for report output
    
    // disable default, copy constructors
    AltChebAlg();
    AltChebAlg(AltChebAlg<Scalar> const &);
    
    
  };

  /** policy class for creation of ChebAlg - build
      method creates ChebAlg with these attributes:

      nrtol    = normal residual (LS gradient) tolerance for convergence
      maxcount = max number of iterations permitted

      Default values set to cause immediate return from ChebAlg::run.

      Other attributes are arguments of build method.

      Conforms to specifications described in TRGNAlg docs.

      Usage: use as policy type, i.e. class name as template parameter
      to TRGNAlg constructor. After construction of TRGNAlg, which IS
      a ChebPolicy by inheritance, call the assign method on it to set
      rtol, nrtol, and maxcount, BEFORE calling TRGNAlg::run() - else
      you will get immediate termination, as intended.
  */
    
    
  template<typename Scalar>
  class AltChebPolicyData {
    
    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;
    
  public:
    
    int maxcount;
    atype gamma;     // inversion level
    atype epsilon;   // error reduction
    atype alpha;     // 'fudge factor'
    atype rhoest;   // estimated spectrum bound
    bool verbose;
    
    AltChebPolicyData(atype _gamma = 0.04,
		      atype _epsilon = 0.01,
		      atype _alpha = 1.001,
		      atype _rhoest = 0.0,
		      int _maxcount = 0,
		      bool _verbose = false)
      : maxcount(_maxcount), gamma(_gamma), epsilon(_epsilon), alpha(_alpha),rhoest(_rhoest),verbose(_verbose) {}
    
    AltChebPolicyData(AltChebPolicyData<Scalar> const & a)
      : maxcount(a.maxcount), gamma(a.gamma), epsilon(a.epsilon), alpha(a.alpha), rhoest(a.rhoest), verbose(a.verbose) {}
      
    ostream & write(ostream & str) const {
      str<<"\n";
      str<<"==============================================\n";
      str<<"AltChebPolicyData: \n";
      str<<"gamma     = "<<gamma<<"\n";
      str<<"epsilon   = "<<epsilon<<"\n";
      str<<"alpha     = "<<alpha<<"\n";
      str<<"rhoest    = "<<rhoest<<"\n";
      str<<"maxcount  = "<<maxcount<<"\n";
      str<<"verbose   = "<<verbose<<"\n";
      str<<"==============================================\n";
      return str;
    }
    
  };
  
  template<typename Scalar> 
  class AltChebPolicy {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  public:
    /**
       build method - see TRGNAlg specs

       @param x - solution vector, initialize to zero on input,
       estimated solution on output 

       @param opeval - Operator evaluation, from which get linear op
       (getDeriv) and RHS (getValue), to set up Gauss-Newton problem

       @param nrnorm - reference to normal residual norm scalar, norm
       of normal residual (least squares gradient) on input, updated
       to estimated solution on output
       
       @param Delta - trust radius (const) passed by value @param str
       - verbose output stream
    */
    AltChebAlg<Scalar> * build(Vector<Scalar> & x, 
			       LinearOp<Scalar> const & A,
			       Vector<Scalar> const & d,
			       atype & nrnorm,
			       ostream & str) const {
      if (verbose)
	return new AltChebAlg<Scalar>(x, A, d, nrnorm,gamma,epsilon,alpha,rhoest,maxcount,str);
      else
	return new AltChebAlg<Scalar>(x, A, d, nrnorm, gamma,epsilon,alpha,rhoest,maxcount,nullstr);
    }

    /** post-construction initialization
	@param _nrtol - normal residual (LS gradient) norm stopping threshhold
	@param _maxcount - max number of permitted iterations
    */
    void assign(atype _gamma, atype _epsilon, atype _alpha, atype _rhoest, int _maxcount, bool _verbose) {
      gamma = _gamma; epsilon=_epsilon; alpha=_alpha; rhoest=_rhoest; maxcount=_maxcount; verbose=_verbose;
    }

    /** parameter table overload */
    void assign(Table const & t) {
      gamma=getValueFromTable<atype>(t,"Cheb_gamma");
      epsilon=getValueFromTable<atype>(t,"Cheb_epsilon");
      alpha=getValueFromTable<atype>(t,"Cheb_alpha");
      rhoest=getValueFromTable<atype>(t,"Cheb_rhoest");
      maxcount=getValueFromTable<int>(t,"Cheb_MaxItn");
      verbose=getValueFromTable<bool>(t,"Cheb_Verbose");
    }

    /** data struct overload */
    void assign(AltChebPolicyData<Scalar> const & s) {
      gamma=s.gamma;
      epsilon=s.epsilon;
      alpha=s.alpha;
      rhoest=s.rhoest;
      maxcount=s.maxcount;
      verbose=s.verbose;
    }

    /** main constructor - acts as default. Default values of
	parameters set to result in immediate return, no
	iteration. Note that policy design requires that default
	construction must be valid, and all run-time instance data be
	initiated post-construction, in this case by the assign
	function, to be called by drivers of user classes (subclassed
	from this and with this as template param).
    */
    AltChebPolicy(atype _gamma = 0.04,
		  atype _epsilon = 0.01,
		  atype _alpha = 1.001,
		  atype _rhoest = 0.0,
		  int _maxcount = 0,
		  bool _verbose = true)
      : gamma(_gamma), epsilon(_epsilon), alpha(_alpha), rhoest(_rhoest),maxcount(_maxcount), verbose(_verbose), nullstr(0){}
      
    AltChebPolicy(AltChebPolicy<Scalar> const & p)
      : gamma(p.gamma),
	epsilon(p.epsilon),
	alpha(p.alpha),
	rhoest(p.rhoest),
	maxcount(p.maxcount),
	verbose(p.verbose),
	nullstr(0) {}
      
  private:
    mutable atype gamma;
    mutable atype epsilon;
    mutable atype alpha;
    mutable atype rhoest;
    mutable int maxcount;
    mutable bool verbose;
    mutable std::ostream nullstr;
  };
}

#endif









