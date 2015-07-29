// cgalg.H
// created by ADP
// last modified 06/17/04

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

#ifndef __RVL_CGALG
#define __RVL_CGALG

#include "alg.hh"
#include "terminator.hh"
#include "linop.hh"

namespace RVLUmin {

  using namespace RVL;
  using namespace RVLAlg;

  /** used to test curvature in CG - handles complex case properly
   */
  template<typename Scalar>
  bool realgt(Scalar left, Scalar right) {
    if (left > right) return true;
    return false;
  }

  template<typename Scalar>
  bool realgt(complex<Scalar> left, complex<Scalar> right) {
    if (real(left) > real(right)) return true;
    return false;
  }

  /** Exception subtype - thrown when needed */
  class CGException: public RVLException {
  public:
    CGException(): RVLException() {
      (*this)<<"Error: RVLUmin::CGAlg::run\n";
    }
    CGException(CGException const & s): RVLException(s) {}
    ~CGException() throw() {}
  };

  /** Single iteration of the Conjugate Gradient method for solution
      of SPD linear systems. Stores const references for linear map,
      represented as RVL::LinearOp, and right-hand side vector,
      represented as RVL::Vector. Updates solution vector and residual
      norm, for which the class stores mutable references. Detection
      of negative curvature (implying failure of SPD condition) throws
      exception. Uses RVL::ProtectedDivision to compute alpha, beta,
      with default zerodivide tolerance. */
  template<typename Scalar>
  class CGStep: public Algorithm {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;
  
  public:

    /**
       Execute a single CG step.  If CG finds a direction of
       negative curvature <p, Ap> < 0, it will NOT take a step
    */
    void run() {

      Scalar alpha, beta; 

      A.applyOp(p,w);
      curv = p.inner(w);

      if (!realgt(curv,ScalarFieldTraits<Scalar>::Zero()) ) {
	CGException e;
	e<<"RVLUmin::CGStep::run: termination\n";
	e<<"  negative curvature (p^TAp) in search direction p = "<<curv<<"\n";
	throw e;
      }

      if (ProtectedDivision<Scalar>(rnormsq,curv,alpha)) {
	CGException e;
	e<<"RVLUmin::CGStep::run: termination\n";
	e<<"  curvature much smaller than mean square residual, yields zerodivide\n";
	e<<"  curvature = "<<curv<<"\n";
	e<<"  mean square residual  = "<<rnormsq<<"\n";
	throw e;
      }

      x.linComb(alpha, p);    
      r.linComb(-alpha, w);
      beta = rnormsq;
      rnormsq = r.normsq();
      
      if (ProtectedDivision<Scalar>(rnormsq,beta,beta)) {
	CGException e;
	e<<"RVLUmin::CGStep::run: termination\n";
	e<<"  previous square residual much smaller than current, yields zerodivide\n";
	e<<"  previous square residual = "<<beta<<"\n";
	e<<"  current square residual  = "<<rnormsq<<"\n";
	throw e;
      }

      p.linComb(ScalarFieldTraits<Scalar>::One(), r, beta);
    }
     
    /** constructor: 

	@param x0 - mutable reference to solution vector (external),
	updated by CGStep::run()

	@param inA - const reference to LinearOp (external) defining
	problem; presumed to be SPD

	@param rhs - const reference to RHS or target vector
	(external)

	@param _rnormsq - mutable reference to residual norm-squared scalar;
	constructor initializes, CGStep::run() updates

    */
    CGStep( Vector<Scalar> & x0, LinearOp<Scalar> const & inA, 
	    Vector<Scalar> const & rhs, atype & _rnormsq)
      : x(x0), A(inA), b(rhs), r(A.getRange()), 
	curv(numeric_limits<Scalar>::max()),
	rnormsq(_rnormsq), 
	p(A.getRange()), 
	w(A.getRange()) {
      CGStep<Scalar>::restart();
    }

 protected:
    
    /** restart CG for the current x.
	This should only be used if x is modified outside of this
	algorithm.
    */
    void restart() {
      A.applyOp(x, w);
      r.copy(b);
      r.linComb(-1.0, w);
      rnormsq = r.normsq();
      p.copy(r);
    }

    Vector<Scalar> & x;         // the current iterate
    const LinearOp<Scalar> & A; // the linear op we want to invert
    Vector<Scalar> const & b;   // the right hand side
    Vector<Scalar> r;           // r = b - A * x
    Scalar curv;                // p' * A * p
    atype & rnormsq;            // a reference to residual norm squared = r' * r
    Vector<Scalar> p;           // the orthogonal part of the residual
    Vector<Scalar> w;           // A * p
  };

  /** implementation of a CG algorithm. Combines CGStep
      with a terminator which displays iteration count and mean square
      residual, and terminates if iteration count exceeds max or
      residual norm tolerance falls below threshhold (default =
      10*sqrt(macheps)). Note that the thresshold is expressed in
      terms of the norm squared, not the norm. Also terminates if step
      (net, from initial estimate of solution) exceeds maxstep
      argument to constructor. In this latter case, the computed step
      is projected onto the ball of radius maxstep centered at the
      initial estimate. This maximum step limit and projection turns
      the algorithm into an approximate trust region subproblem
      solver. The default choice of maxstep is the max Scalar, which
      effectively turns off the trust region feature.

      Note that the solution vector (constructor parameter _x) is not
      altered by the constructor, and represents an externally defined
      initial estimate. In particular, it is not presumed to be the
      zero vector, nor is it set to zero by the constructor.

      For definition of arguments to constructor, see below. 

      Usage: construct CGAlg object, supplying arguments as indicated below. Invoke run() method.

      Typical use case: see <a href="../../testsrc/testcg.cc">functional test source</a>.
  */
  template<typename Scalar>
  class CGAlg: public Algorithm {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  public:

    /** constructor: 

	@param _x - mutable reference to solution vector (external),
	estimated solution on return from CGAlg::run().

	@param _inA - const reference to LinearOp (external) defining
	problem; presumed to be SPD

	@param _rhs - const reference to RHS or target vector
	(external)

	@param _rnormsq - mutable reference to residual norm-squared scalar
	(external), initialized to norm-squared of RHS on construction, norm-squared
	of estimated residual at solution on return from
	CGAlg::run()

	@param _tol - stopping threshold for residual norm-squared, default
	value = 100.0*macheps

	@param _maxcount - max number of iterations permitted, default
	value = 10

	@param _maxstep - max permitted step length (trust radius),
	default value = max absval scalar (which makes the trust
	region feature inactive)

	@param _str - output stream
    */
    CGAlg(RVL::Vector<Scalar> & _x, 
		LinearOp<Scalar> const & _inA, 
		Vector<Scalar> const & _rhs, 
		atype & _rnormsq,
		atype _tol = 100.0*numeric_limits<atype>::epsilon(),
		int _maxcount = 10,
		atype _maxstep = numeric_limits<atype>::max(),
		ostream & _str = cout)
      : x(_x),
	resname("MS residual"), 
	inA(_inA), 
	rhs(_rhs), 
	tol(_tol), 
	maxcount(_maxcount), 
	maxstep(_maxstep), 
	str(_str), 
	rnormsq(_rnormsq),
	step(x,inA,rhs,rnormsq),
	it(0) {}

    int getCount() { return it; }
    void run() { 
      try {
	// terminator for CG iteration
	CountingThresholdIterationTable<Scalar> stop1(maxcount,rnormsq,tol,resname,str);
	// terminator for Trust Region test and projection
	BallProjTerminator<Scalar> stop2(x,maxstep,str);
	// terminate if either
	OrTerminator stop(stop1,stop2);
	// loop
	LoopAlg doit(step,stop);
	doit.run();
	
	// must recompute residual if scaling occured 
	if (stop2.query()) {
	  str<<"CGAlg::run: scale step to trust region boundary\n";
	  Vector<Scalar> temp(inA.getDomain());
	  inA.applyOp(x,temp);
	  temp.linComb(-1.0,rhs);
	  rnormsq=temp.normsq();
	}

	it=stop1.getCount();
      }
      catch (CGException & e) {
	throw e;
      }
      catch (RVLException & e) {
	e<<"Error: CGAlg::run\n";
	throw e;
      }
    }

  private:

    LinearOp<Scalar> const & inA;  // operator - should be SPD
    Vector<Scalar> & x;            // state - solution vector
    Vector<Scalar> const & rhs;    // reference to rhs
    atype tol;                     // tolerance for norm squared of residual
    atype maxstep;                 // upper bound for net step x-x0
    int maxcount;                  // upper bound for iteration count
    string resname;                // string to print in output table
    ostream & str;                 // stream for report output
    int it;                        // private iteration count
    atype & rnormsq;               // residual norm squared
    CGStep<Scalar> step;           // alg expressing single CG step

  };

}

#endif
