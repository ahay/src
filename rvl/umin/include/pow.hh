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

#ifndef __RVL_POWALG
#define __RVL_POWALG

#include "alg.hh"
#include "linop.hh"
#include "ioterm.hh";

namespace RVLUmin {

  using namespace RVL;
  using namespace RVLAlg;

  /** This Algorithm does a single iteration of the Power Method for
      estimating the largest singular value of a linear operator.  It
      also estimates the relative residual |lam -rq|/|lam|, in which
      lam is the eigenvalue of the normal operator which the Rayleigh
      quotient rq approximates.
  */
  template <class Scalar>
  class PowerStep: public Algorithm {

    typedef typename ScalarFieldTraits<Scalar>::AbsType AScalar;

  private:

    Vector<Scalar> & x;         // initial guess at eigenvector on construction, current iterate else
    const LinearOp<Scalar> & A; // the linear op 
    Vector<Scalar> p;           // workspace for A * x
    Vector<Scalar> w;           // workspace for A^T * p
    AScalar & sig;              // estimated largest singular value
    AScalar & err;              // residual norm
    AScalar rq;                 // Rayleigh quotient
    AScalar sx;                 // workspace for 1/x.norm

    PowerStep();
    PowerStep(PowerStep< Vector<Scalar> > const &);

  public:
  
    /** Constructor parameters:
	@param _x = estimated singular vector - mutable reference, updated on return
	@param _A = linear operator for which singular pair is sought - const reference
	@param _sig = estimated max singular value - mutable reference, updated on return
	@param _err = estimate of relative error in singular value - mutable reference, updated on return
    */
    PowerStep( RVL::Vector<Scalar> & _x, const LinearOp<Scalar> & _A,
	     AScalar & _sig, AScalar & _err )
      : x(_x), 
	A(_A), 
	p(A.getDomain()), 
	w(A.getRange()),
	sig(_sig),
	err(_err),
	sx(ScalarFieldTraits<AScalar>::One())
    {	  
      AScalar one = ScalarFieldTraits<AScalar>::One();
      AScalar sx = one;
      if (ProtectedDivision<AScalar>(one,x.norm(),sx)) {
	RVLException e;
	e<<"Error: RVLUmin::PowerStep::constructor\n";
	e<<"zerodivide\n";
	throw e;
      }
      x.scale(sx);
    }

    /**
       Execute a single power step.
    */
    void run() {

      try {
	AScalar one = ScalarFieldTraits<AScalar>::One();
	A.applyOp(x,w);
	rq = w.normsq();
	if (rq < ScalarFieldTraits<AScalar>::Zero()) {
	  RVLException e;
	  e<<"Error: PowerMethod::run\n";
	  e<<"computed apparently negative Rayleigh quotient\n";
	  throw e;
	}
	sig=sqrt(rq);
	if (ProtectedDivision<AScalar>(one,rq,err)) {
	  sig=ScalarFieldTraits<AScalar>::Zero();
	  err=ScalarFieldTraits<AScalar>::One();
	}
	else {
	  A.applyAdjOp(w,p);
	  // p <- p-rq*x
	  p.linComb(-rq,x);
	  // err=|p-rq*x|/|rq| \simeq |\lambda-rq|/|rq| \simeq |\lambda-rq|/|\lambda|
	  err*=p.norm();
	  // p <- p+rq*x - restores p = A^TAx, avoids another Vector workspace
	  p.linComb(rq,x);
	  // compute normalizing factor
	  if (ProtectedDivision<AScalar>(one,p.norm(),sx)) {
	    RVLException e;
	    e<<"Error: RVLUmin::PowerStep::run\n";
	    e<<"zerodivide\n";
	    throw e;
	  }
	  // x = A^TAx / |A^TAx|
	  x.scale(sx,p);
	}
      }

      catch (RVLException & e) {
	e<<"\ncalled from PowerStep::run\n";
	throw e;
      }
    }
  };

  /** Power method for finding largest singular value of a linear operator. 

      Algorithm: stopping criterion based on an estimate of the
      relative residual \f$|\lambda - \rho|/|\lambda|\f$, in which
      \f$\lambda\f$ is the largest eigenvalue of the normal operator
      and \f$\rho\f$ is the current Rayleigh quotient. Since
      \f$\lambda\f$ is not available, this numerator is approximated
      by the residual norm, the denominator by \f$\rho\f$: thus the
      error estimator is \f$|A^tAx - \rho x|/|\rho|\f$.

      Parameters: The only two characteristic parameters required by
      this algorithm are the max number of steps permitted, and the
      stopping tolerance for the relative residual test. Reasonable
      values for these depend on the convergence rate, which is a
      function of the distribution of singular values. An isolated
      largest singular value (or group of singular values near the
      max) produces fastest convergence, with an elementary estimate
      showing dependence of the rate on the gap between largest and
      next-largest singular values (eg. Parlett, <i>The Symmetric
      Eigenvalue Problem</i>, SIAM, 1998). Thus, no general guidance
      may be given for selection of these parameters.

      Other parameters supplied to the constructor are mutable
      references for the singular value and vector workspace, and a
      const reference for the operator. See constructor docs.

      Typical use case: see <a
      href="../../testsrc/testpow.cc">functional test source</a>.

   */
  template<typename Scalar>
  class PowerMethod: public Algorithm {

    typedef typename ScalarFieldTraits<Scalar>::AbsType AScalar;

  private:

    AScalar err;                        // workspace for rel residual calc
    PowerStep<Scalar> step;             // step alg, statically allocated
    std::vector<std::string> names;     // headings for report cols
    std::vector<AScalar *> nums;        // pointers to err, sig
    std::vector<AScalar> tols;          // convergence tol for err, zero for sig (inactive)
    int _nstep;                         // max permitted steps
    int it;                             // convenience store of steps taken, for tests
    ostream & _str;                     // output stream

  public:
    
    /** Constructor parameters:
	@param x = initial singular vector estimate - mutable reference, updated on return
	@param sig = initial singular value estimate - mutable reference, updated on return
	@param A = linear operator for which singular pair is sought - const reference
	@param nstep = number of power method iterations permitted. Typical value = 50.
	@param tol = stopping threshold for error estimator. Typical value = 0.001
	@param str = output stream on which to report progress of iteration
    */
    PowerMethod(Vector<Scalar> & x,
		AScalar & sig,
		LinearOp<Scalar> const & A,
		int nstep, AScalar tol, ostream & str) 
      : err(ScalarFieldTraits<AScalar>::One()),
	step(x,A,sig,err),_nstep(nstep),_str(str),
	names(2),nums(2),tols(2), it(0) {
      names[0]="Relative Residual"; nums[0]=&err; tols[0]=tol;
      names[1]="Est Singular Value"; nums[1]=&sig; tols[1]=ScalarFieldTraits<AScalar>::Zero();
    }
    
    void run() {
      try {
	VectorCountingThresholdIterationTable<Scalar> term(_nstep,names, nums, tols,_str);
	term.init();
	LoopAlg l(step,term);
	l.run();
	it=term.getCount();
      }
      catch (RVLException & e) {
	e<<"\ncalled from PowerMethod::run\n";
	throw e;
      }
    }
    
    int getCount() { return it; }
  };

}

#endif
