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

  /** Single step of Chebyshev iteration for the normal
      equations.
      
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
  class ChebStep : public Algorithm {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  public:

    ChebStep(LinearOp<Scalar> const & _A,
	     Vector<Scalar> & _x,
	     Vector<Scalar> const & _b,
	     atype & _rnorm, 
	     atype & _nrnorm,
         std::vector<atype> &_coeff,
         int & _kc,
         int & _nrt,
         atype _gamma = 0.04,     // inversion level
         atype _epsilon = 0.001,  // error reduction
         atype _alpha = 1.1,      // 'fudge factor'
         atype _lbd_est = 0.0,    // input operator spectrum in case it is know prior
         ostream & _str = cout)
      : A(_A), x(_x), b(_b), rnorm(_rnorm), nrnorm(_nrnorm), coeff(_coeff), gamma(_gamma), epsilon(_epsilon), alpha(_alpha), lbd_est(_lbd_est),kc(_kc), nrt(_nrt), str(_str), dx(A.getDomain()), r(A.getDomain()), ndx(A.getDomain()), tmp_vector(A.getRange()), x_init(_x),dx_init(A.getDomain()),r_init(A.getDomain()),ndx_init(A.getDomain()){
        
      atype one = ScalarFieldTraits<atype>::One();
      atype tmp = one;
      nrt = 0;
      // NOTE: set up dx_init, r_init, ndx_int
      // The real step of initialization is in the run() method
      rnorm = b.norm();
      A.applyAdjOp(b,r_init);
      dx_init.zero();
      nrnorm=r_init.norm();
      ndx_init.zero();
      A.applyOp(r_init,tmp_vector);
      A.applyAdjOp(tmp_vector,ndx_init);
      //nrnorm=ndx_init.norm();
      
      tmp = abs(r_init.inner(ndx_init));
      if (ProtectedDivision<atype>(tmp,r_init.normsq(),RQ)) {
          RVLException e;
          e<<"Error: ChebStep::ChebStep() from ProtectedDivision: RQ\n";
          throw e;
      }
      if(RQ > lbd_est){
          lbd_est = alpha * RQ;
          nrt = nrt + 1;
      }

      ProtectedDivision<atype>(2*one,(one+gamma)*lbd_est,s);
    }
      
    /**
       Run a single step of the Chebyshev iteration for the normal equations
    */
    void run() {
      try {
          atype one = ScalarFieldTraits<atype>::One();
          if (kc == 0){
              x.copy(x_init);
              r.copy(r_init);   
              dx.copy(dx_init);
              ndx.copy(ndx_init);
              str<<"Estimated spectrum bound at iter ["<<nrt<<"] = "<< lbd_est << endl;
              str << "NOTE: The  " << nrt << "-th restart\n";
          }
      
          atype tmp;
          dx.linComb(s*coeff[kc+1],r,coeff[kc+1]-one); // this = a*x+b*this
          
          x.linComb(one,dx);
          A.applyOp(dx,tmp_vector);
          A.applyAdjOp(tmp_vector,ndx);
          A.applyOp(x,tmp_vector);
          tmp_vector.linComb(-1*one,b);
          rnorm = tmp_vector.norm();
          
          r.linComb(-1*one,ndx);
          nrnorm=r.norm();
          //nrnorm = ndx.norm();
          tmp = abs(dx.inner(ndx));
          if (ProtectedDivision<atype>(tmp,dx.normsq(),RQ)) {
              RVLException e;
              e<<"Error: ChebStep::run() from ProtectedDivision: RQ\n";
              throw e;
          } 
          if(RQ > lbd_est){
              lbd_est = alpha * RQ;
              ProtectedDivision<atype>(2*one,(one+gamma)*lbd_est,s);
              kc = -1;
              nrt = nrt + 1;
          }
          kc = kc + 1;
      }
      catch (RVLException & e) {
	e<<"\ncalled from ChebStep::run()\n";
	throw e;
      }
     
    }
    
    atype getSpectrumBound() const { return lbd_est; }

    ~ChebStep() {}

  private:

    // references to external objects
    LinearOp<Scalar> const & A;
    Vector<Scalar> & x;
    Vector<Scalar> const & b;
    atype & rnorm;
    atype & nrnorm;
    
    // added for Chebyshev
    std::vector<atype> &coeff;
 
    atype gamma; 		   // inversion level
    atype epsilon; 		   // error reduction factor
    atype alpha;		   // 'fudge factor'
    atype beta;
    atype RQ;
    atype lbd_est;         // estimated spectrum bound
    atype s;
    int &kc;
    int &nrt;
    ostream & str;
    // end of added
    
    // need four work vectors and one scalar as persistent object data
    Vector<Scalar> dx;
    Vector<Scalar> r;
    Vector<Scalar> ndx;
    Vector<Scalar> tmp_vector;
      
    // need vectors to store initial values for restarting purpose
    Vector<Scalar> x_init;
    Vector<Scalar> dx_init;
    Vector<Scalar> r_init;
    Vector<Scalar> ndx_init;
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
  class ChebAlg: public Algorithm {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  public:

    /** Constructor

	@param _x - mutable reference to solution vector (external),
	initialized to zero vector on construction, estimated solution
	on return from ChebAlg::run().

	@param _inA - const reference to LinearOp (external) defining
	problem

	@param _rhs - const reference to RHS or target vector
	(external)

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
    ChebAlg(RVL::Vector<Scalar> & _x, 
	    LinearOp<Scalar> const & _inA, 
	    Vector<Scalar> const & _rhs, 
	    atype & _rnorm,
	    atype & _nrnorm,
	    atype _gamma = 0.04,     // inversion level
	    atype _epsilon = 0.001,  // error reduction
	    atype _alpha = 1.1,      // 'fudge factor'
            atype _lbd_est=0.0,
	    int _maxcount = 10,      // upper bound of iterations
	    ostream & _str = cout)  
      : inA(_inA), x(_x), rhs(_rhs), rnorm(_rnorm), nrnorm(_nrnorm), gamma(_gamma),epsilon(_epsilon),alpha(_alpha),lbd_est(_lbd_est),maxcount(_maxcount), kmax(_maxcount), kc(0), ktot(0), str(_str), step(inA,x,rhs,rnorm,nrnorm,coeff,kc,nrt,gamma,epsilon,alpha,lbd_est,str)
	// NOTE: reference to coeff has been passed to ChebStep object step, however
	// coeff has not been initialized.
    { x.zero();
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
    str << "NOTE: computed number of iterations needed:  " << kmax << "\n";
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
	  //str<<"coeff["<<i<<"] = "<< coeff[i] << endl;
	}
      }
  

    ~ChebAlg() {}

    void run() { 
      
      // terminator for Cheb iteration
      vector<string> names(2);
      vector<atype *> nums(2);
      vector<atype> tols(2);
      atype rnorm0=rnorm;
      atype nrnorm0=nrnorm;
      names[0]="Residual Norm"; nums[0]=&rnorm; tols[0]=ScalarFieldTraits<atype>::Zero();
      names[1]="Normal Residual Norm"; nums[1]=&nrnorm; tols[1]=tols[0]=ScalarFieldTraits<atype>::Zero();
      str<<"========================== BEGIN Cheb =========================\n";
      VectorCountingThresholdIterationTable<atype> stop1(maxcount,names,nums,tols,str);
      MaxTerminator<int> stop2(kc,kmax-1);
      //MaxTerminator<int> stop3(nrt,maxrt);
      stop1.init();
      // terminate if either
      OrTerminator stop(stop1,stop2);
      // loop
      LoopAlg doit(step,stop);
      doit.run();
      ktot = stop1.getCount();

      str<<"=========================== END Cheb ==========================\n";
        // display results
        str<<"\n ******* summary ********  "<<endl;
        str<<"initial residual norm      = "<<rnorm0<<endl;
        str<<"residual norm              = "<<rnorm<<endl;
        str<<"residual redn              = "<<rnorm/rnorm0<<endl;
        str<<"initial gradient norm      = "<<nrnorm0<<endl;
        str<<"gradient norm              = "<<nrnorm<<endl;
        str<<"gradient redn              = "<<nrnorm/nrnorm0<<endl;
    }

    int getCount() const { return kc; }
    int getTotalCount() const { return ktot; }
    int getRestartCount() const { return nrt+1; }
    atype getSpectrumBound() const { return step.getSpectrumBound(); }


  private:

    LinearOp<Scalar> const & inA;  // operator
    Vector<Scalar> & x;            // state - solution vector
    Vector<Scalar> const & rhs;    // reference to rhs
    atype & rnorm;                 // residual norm
    atype & nrnorm;                // gradient norm

    // added for Chebyshev
    std::vector<atype> coeff;
    atype gamma; 		   // inversion level
    atype epsilon; 		   // error reduction factor
    atype alpha;		   // 'fudge factor'
    atype lbd_est;         // spectrum bound of op
    int maxcount;                  // upper bound for iteration count
    int nrt;                       // upper bound for restarting
    // end of added 

    /** Necessary number of iterations for a given error reduction factor.
	Special for Chebyshev iteration*/
    int kmax;
    int kc;
    int ktot;
    ostream & str;                 // stream for report output
    ChebStep<Scalar> step;         // single step of Cheb

    // disable default, copy constructors
    ChebAlg();
    ChebAlg(ChebAlg<Scalar> const &);

  };

  /** policy class for creation of ChebAlg - build
      method creates ChebAlg with these attributes:

      rtol     = residual threshhold for convergence
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
  class ChebPolicyData {
    
    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;
    
  public:
    
    int maxcount;
    atype gamma;     // inversion level
    atype epsilon;   // error reduction
    atype alpha;     // 'fudge factor'
    atype lbd_est;   // estimated spectrum bound
    bool verbose;
    
    ChebPolicyData(atype _gamma = 0.04,
		   atype _epsilon = 0.01,
		   atype _alpha = 1.001,
           atype _lbd_est = 0.0,
		   int _maxcount = 0,
           bool _verbose = false)
      : maxcount(_maxcount), gamma(_gamma), epsilon(_epsilon), alpha(_alpha),lbd_est(_lbd_est),verbose(_verbose) {}
    
    ChebPolicyData(ChebPolicyData<Scalar> const & a)
      : maxcount(a.maxcount), gamma(a.gamma), epsilon(a.epsilon), alpha(a.alpha), lbd_est(a.lbd_est), verbose(a.verbose) {}
      
    ostream & write(ostream & str) const {
          str<<"\n";
          str<<"==============================================\n";
          str<<"ChebPolicyData: \n";
          str<<"gamma     = "<<gamma<<"\n";
          str<<"epsilon   = "<<epsilon<<"\n";
          str<<"alpha     = "<<alpha<<"\n";
          str<<"lbd_est   = "<<lbd_est<<"\n";
          str<<"maxcount  = "<<maxcount<<"\n";
          str<<"verbose   = "<<verbose<<"\n";
          str<<"==============================================\n";
          return str;
      }
    
  };
  
  template<typename Scalar> 
  class ChebPolicy {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  public:
    /**
       build method - see TRGNAlg specs

       @param x - solution vector, initialize to zero on input,
       estimated solution on output 

       @param opeval - Operator evaluation, from which get linear op
       (getDeriv) and RHS (getValue), to set up Gauss-Newton problem

       @param rnorm - reference to residual norm scalar, norm of RHS
       on input, of residual on output

       @param nrnorm - reference to normal residual norm scalar, norm
       of normal residual (least squares gradient) on input, updated
       to estimated solution on output
       
       @param Delta - trust radius (const) passed by value @param str
       - verbose output stream
    */
    ChebAlg<Scalar> * build(Vector<Scalar> & x, 
			    LinearOp<Scalar> const & A,
			    Vector<Scalar> const & d,
			    atype & rnorm,
			    atype & nrnorm,
			    ostream & str) const {
        if (verbose)
      return new ChebAlg<Scalar>(x, A, d,rnorm,nrnorm,gamma,epsilon,alpha,lbd_est,maxcount,str);
        else
      return new ChebAlg<Scalar>(x, A, d,rnorm,nrnorm,gamma,epsilon,alpha,lbd_est,maxcount,nullstr);
    }

    /** post-construction initialization
	@param _rtol - residual norm stopping threshhold
	@param _nrtol - normal residual (LS gradient) norm stopping threshhold
	@param _maxcount - max number of permitted iterations
    */
    void assign(atype _rtol, atype _nrtol, atype _gamma, atype _epsilon, atype _alpha, atype _lbd_est, int _maxcount, bool _verbose) {
        gamma = _gamma; epsilon=_epsilon; alpha=_alpha; lbd_est=_lbd_est; maxcount=_maxcount; verbose=_verbose;
    }

    /** parameter table overload */
    void assign(Table const & t) {
      gamma=getValueFromTable<atype>(t,"Cheb_gamma");
      epsilon=getValueFromTable<atype>(t,"Cheb_epsilon");
      alpha=getValueFromTable<atype>(t,"Cheb_alpha");
      lbd_est=getValueFromTable<atype>(t,"Cheb_lbd_est");
      maxcount=getValueFromTable<int>(t,"Cheb_MaxItn");
      verbose=getValueFromTable<bool>(t,"Cheb_Verbose");
    }

    /** data struct overload */
    void assign(ChebPolicyData<Scalar> const & s) {
      gamma=s.gamma;
      epsilon=s.epsilon;
      alpha=s.alpha;
      lbd_est=s.lbd_est;
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
    ChebPolicy(atype _gamma = 0.04,
               atype _epsilon = 0.01,
               atype _alpha = 1.001,
               atype _lbd_est = 0.0,
	       int _maxcount = 0,
            bool _verbose = true)
      : gamma(_gamma), epsilon(_epsilon), alpha(_alpha), lbd_est(_lbd_est),maxcount(_maxcount), verbose(_verbose), nullstr(0){}
      
    ChebPolicy(ChebPolicy<Scalar> const & p)
      : gamma(p.gamma),
    epsilon(p.epsilon),
    alpha(p.alpha),
    lbd_est(p.lbd_est),
    maxcount(p.maxcount),
    verbose(p.verbose),
    nullstr(0) {}
      
  private:
    mutable atype gamma;
    mutable atype epsilon;
    mutable atype alpha;
    mutable atype lbd_est;
    mutable int maxcount;
    mutable bool verbose;
    mutable std::ostream nullstr;
  };
}

#endif









