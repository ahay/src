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

#ifndef __RVLALG_UMIN_LBFGSALG_H
#define __RVLALG_UMIN_LBFGSALG_H

#include "newtonalg.hh"
#include "umintable.hh"
#include "terminator.hh"

namespace RVLUmin{

  using namespace RVLAlg;
  using RVL::Vector;
  using RVL::LinearOp;
  using RVL::Functional;
  using RVL::FunctionalEvaluation;

  /** {\bf LMBFGSOp} implements the limited memory BFGS approximation
      to the inverse Hessian of a twice-differentiable function.  This
      approximation uses local changes to the gradient to gradually
      build up an estimate of the Hessian for use in nonlinear
      optimization problems.
    
      For details of the algorithm, see the paper
    
      "Updating Quasi-Newton Matrices with Limited Storage" by Jorge Nocedal,
      Math. of Computation, Vol. 35, no. 151, p.p. 773--782.
    
      Note that the operator is an approximation to the {\it inverse} Hessian,
      so the the {\bf apply} method computes the quasi-Newton step.
    
      The BFGS approximation is based on successive rank-one perturbations to
      an initial approximation; these rank-one perturbations can be represented
      as outer products.  This class allows the user to provide a symmetric
      positive definite operator to define an alternate inner product, which
      then changes the definition of the outer product.  In the context of an
      optimization problem, this is equivalent to implementing the algorithm
      in the alternate inner product.  */

  template<class Scalar>
  class LBFGSOp : public LinearOp<Scalar> {  
  
  private:

    const Space<Scalar> & sp;
    const CartesianPowerSpace<Scalar> psp;
    ScaleOpFwd<Scalar> H0;
    int CurNum,  // indicates the last vector allocated
      MaxNum,    // Maxnum == psp.getSize() == # of vectors
      CurAllocated; // if CurAllocated < MaxNum, then empty slots
    Vector<Scalar> Yvec;
    Components<Scalar> Y;
    Vector<Scalar> Svec;
    Components<Scalar> S;
    std::vector<Scalar> rho;

    // default and copy constructors---disabled.
    LBFGSOp();
  
  public:

    /** Usual constructor. Needs a multiple of the identity to use for the
	initial inverse Hessian approximation, the maximum number of
	updates, and, optionally, an operator to change the inner product. */
    LBFGSOp(const Space<Scalar> & _sp,
	    Scalar ihs,
	    int maxnum)
      : sp(_sp), psp(maxnum, sp), H0(sp,ihs),CurNum(0),MaxNum(maxnum),
	CurAllocated(0),Yvec(psp), Y(Yvec),Svec(psp),S(Svec),rho(MaxNum) {
    }

    /** Copy constructor */
    LBFGSOp(const LBFGSOp<Scalar> & op) 
      : sp(op.sp), psp(op.MaxNum, sp), H0(op.H0),CurNum(op.CurNum),MaxNum(op.MaxNum),
	CurAllocated(op.CurAllocated),Yvec(psp), Y(Yvec),Svec(psp),S(Svec),rho(op.rho) 
    {
      Yvec.copy(op.Yvec);
      Svec.copy(op.Svec);
    }
  
    // Destructor.
    ~LBFGSOp() {
    }

    const Space<Scalar> & getDomain() const { return sp; }
    const Space<Scalar> & getRange() const { return sp; }

    // Access to the scale, which is dynamically changed by the operator.
    Scalar getScale() {
      try {
	return H0.getScale();
      } 
      catch (RVLException & e) {
	e<<"\ncalled from LBFGSOp::getScale\n";
	throw e;
      }
    }
  
    /** update requires current and next $x$ and current and next gradient. */
    bool update(const Vector<Scalar> & x,
		const Vector<Scalar> & xnext,
		const Vector<Scalar> & g,
		const Vector<Scalar> & gnext) {
      try {
	if (MaxNum==0)
	  return true;
	if (!x.inSpace(sp) || !xnext.inSpace(sp) ||
	    !g.inSpace(sp) || !gnext.inSpace(sp)) {
	  RVLException e; e<<"Error:LBFGSOp::update\n";
	  e<<"input vector(s) not in domain\n";
	  e<<"input vectors:\n";
	  x.write(e); xnext.write(e);
	  g.write(e); gnext.write(e);
	  e<<"domain:\n";
	  sp.write(e);
	  throw e;
	}
      
	if( CurAllocated < MaxNum ) {
	  CurAllocated++;
	}
	S[CurNum].copy(xnext);
	S[CurNum].linComb(-1.0,x);
	Y[CurNum].copy(gnext);
	Y[CurNum].linComb(-1.0,g);
	
	if (ProtectedDivision<Scalar>
	    (1.0,S[CurNum].inner(Y[CurNum]),rho[CurNum])) return false;
	Scalar tmp;
	if (ProtectedDivision<Scalar>
	    (1.0,rho[CurNum]*(Y[CurNum].normsq()),tmp)) return false;
	H0.setScale(tmp);
	CurNum++;
	if(CurNum == MaxNum) CurNum = 0;
	return true;
      }
      catch (RVLException & e) {
	e<<"\ncalled from LBFGSOp::update\n";
	throw e;
      }
    }

    // reset sets the operator to the initial inverse Hessian approximation.
    void reset() { CurNum = 0; CurAllocated = 0; }

    LinearOp<Scalar> * clone() const {
      return new LBFGSOp<Scalar>(*this);
    }

    // apply computes the image of the operator on x, giving y.
    void apply(const Vector<Scalar> & x,
	       Vector<Scalar> & y) const {
      try {
	if (!x.inSpace(sp)) { 
	  RVLException e; e<<"Error:LBFGSOp::apply\n";
	  e<<"input vector not in domain\n";
	  throw e;
	}
	if (!y.inSpace(sp)) {
	  RVLException e; e<<"Error:LBFGSOp::apply\n";
	  e<<"output vector not in range\n";
	  throw e;
	}
      
	Scalar xscale(H0.getScale());

	// handle the special case of no updates
	if (CurAllocated==0) {
	  y.scale(xscale,x);
	  return;
	}

	// general case---the initial approximation has been updated
	std::vector<Scalar> alpha(CurAllocated);
	int i;
	y.copy(x);
	for (i=CurNum-1;i>=0;--i) {
	  alpha[i] = rho[i]*(S[i].inner(y));
	  y.linComb(-alpha[i],Y[i]);
	}
	if( CurAllocated > CurNum )
	  for (i=CurAllocated-1;i>=CurNum;--i) {
	    alpha[i] = rho[i]*(S[i].inner(y));
	    y.linComb(-alpha[i],Y[i]);
	  }

	y.scale(xscale);

	if( CurAllocated > CurNum )
	  for (i=CurNum;i<CurAllocated;++i) {
	    Scalar beta = rho[i]*(Y[i].inner(y));
	    y.linComb(alpha[i]-beta,S[i]);
	  }

	for (i=0;i<CurNum;++i) {
	  Scalar beta = rho[i]*(Y[i].inner(y));
	  y.linComb(alpha[i]-beta,S[i]);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from LBFGSOp::apply\n";
	throw e;
      }
    }

    //  applyAdj computes the image of the adjoint on y, giving x.
    void applyAdj(const Vector<Scalar> & y,
		  Vector<Scalar> & x) const {
      try {
	// operator is by definition symmetric
	apply(y,x);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LBFGSOp::applyAdj\n";
	throw e;
      }
    }

    /// write methods to print out useful information about the object.
    ostream & write(ostream & str) const {
      str<<"LBFGS operator: current rank = "<<CurNum<<" scale = "<<H0.getScale()<<"\n";
      return str;
    }
  };

  /**  This algorithm performs a quasi-newton method for minimizing
       a continuous function.  It uses a limited--memory BFGS approximation
       to the inverse hessian, as formed in the LBFGSOp

       A parameter file may be used to select non--default settings
       for the many parameters.  The name of this file is passed
       to the constructor.
  */

  template<class Scalar>
  class LBFGSStep : public NewtonStep<Scalar> {

  private:

    Vector<Scalar> x0;
    Scalar firststep;
    ostream & str;

    UMinTable<Scalar> const & holder;

    LineSearchAlg<Scalar> const & ls;
    LBFGSOp<Scalar> H;

    LBFGSStep();

  protected:

    // Fill in the vector with the search direction
    virtual bool calcDir(Vector<Scalar> & dir) {
      try {
	// why isn't H defined to be -H?
	const Vector<Scalar> & grad 
	  = ConOptStep<Scalar>::getFunctionalEvaluation().getGradient();
	H.applyOp(grad,dir);
	dir.negate();
	return true;
      } catch(RVLException & e) {
	e << "called from LBFGSStep::calcDir()\n";
	throw e;
      }
      return false;
    }

    /** Find the step length in this search direction,
	and updates x.
	Default implementation always sets the step to 1.
	
	Returning a non-positive step length is 
	an algorithmic failure, and will cause the run()
	method to return false.
    */
    //    virtual bool calcStep(Vector<Scalar> & dir, Scalar & step) {
    virtual bool calcStep(Vector<Scalar> & dir) {
      try {
	//Line Search
	bool tried_steepest_descent = false;

	FunctionalEvaluation<Scalar> & fx =
	  ConOptStep<Scalar>::getFunctionalEvaluation();
	//	cerr<<"fx.getGradient in LBFGSStep\n";
	const Vector<Scalar> & grad = fx.getGradient();
	if (dir.inner(grad) > 0) {
	  str<<"LBFGSStep: encountered nondecrease direction\n";
	  str<<"reset search direction to steepest descent, restart\n";
	  str<<"inverse Hessian computation\n";
	  dir.copy(grad);
          dir.negate();
	  H.reset();
	  tried_steepest_descent = true;
	}

	// set up and run line search
	bool lsr = false;
	bool hur = false;
	ls.initialize(fx,dir,firststep);
	lsr = ls.run();
        if (lsr) { 
	  hur = H.update(ls->getBasePoint(), 
			 fx.getPoint(),
			 ls->getBaseGradient(), 
			 fx.getGradient());
	}
	if (!lsr) {
	  str<<"LBFGSStep: line search failed\n";
	  if (!tried_steepest_descent) {
	    str<<"LBFGSStep: attempting steepest descent restart\n";
	    dir.copy(grad);
	    dir.negate();
	    H.reset();
	    ls.initialize(fx,dir,firststep);
	    lsr = ls.run();
	    if (lsr) { 
	      hur = H.update(ls->getBasePoint(), 
			     fx.getPoint(),
			     ls->getBaseGradient(), 
			     fx.getGradient());
	    }
	  }
	}
	return (lsr && hur);
      } catch(RVLException & e) {
	e << "called from LBFGSStep::calcStep()\n";
	throw e;
      }
    }

  public:

    LBFGSStep(LineSearchAlg<Scalar> const & _ls,
	      FunctionalEvaluation<Scalar> & _fx,
	      Scalar _FirstStep = 1.0,
	      Scalar InvHessianScale = 1.0,
	      int MaxUpdates = 5,
	      ostream & _str = cout)
      : NewtonStep<Scalar>(_fx.getPoint(), _fx), 
	x0(_fx.getPoint()), 
	ls(_ls),
	firststep(_FirstStep),
	H(_fx.getDomain(), 
	  InvHessianScale, 
	  MaxUpdates),
	str(_str) {}
  
    ~LBFGSStep() {} 

  };

  template<class Scalar>
  class UMinLBFGS : public StateAlg<Vector<Scalar> > {

  private:

    LineSearchAlg<Scalar> & ls;
    FunctionalEvaluation<Scalar> fx;

    // LBFGSStep parameters
    Scalar firststep;
    Scalar invhessianscale;
    int maxupdates;

    // Terminator parameters
    int maxitn;
    Scalar tol;
    int displevel;
    bool dispflag;
    ostream & str;
    
    UMinLBFGS();

  public:

    UMinLBFGS(Functional<Scalar> & f,
	      Vector<Scalar> & x,
	      LineSearchAlg<Scalar> & _ls,
	      int _maxitn = 20,
	      Scalar _tol = 0.01,
	      bool _dispflag = false,
	      ostream & _str = cout,
	      Scalar _firststep = 1.0,
	      Scalar _invhessianscale = 1.0,
	      int _maxupdates = 5
	      )
      : ls(_ls), fx(f,x),
	firststep(_firststep),
	invhessianscale(_invhessianscale),
	maxupdates(_maxupdates),
	maxitn(_maxitn),
	tol(_tol),
	dispflag(_dispflag),
        str(_str) {}

    // preferred constructor
    /* deprecated 29.04.07 */
    UMinLBFGS(Functional<Scalar> & f,
	      Vector<Scalar> & x,
	      LineSearchAlg<Scalar> & _ls,
	      UMinTable<Scalar> & holder, 
	      ostream & _str=cout) 
      : ls(_ls), 
	fx(f,x),
	dispflag(false), 
	str(_str) {
      if ((holder.getValue("FirstStep",firststep)) ||
	  (holder.getValue("MaxItn",maxitn)) ||
	  (holder.getValue("InvHessianScale",invhessianscale)) ||
	  (holder.getValue("MaxUpdates",maxupdates)) ||
	  (holder.getValue("GradTol",tol)) ||
	  (holder.getValue("DispFlag",displevel))) {
	RVLException e;
	e<<"Error: UMinLBFGS constructor\n";
	e<<"holder corrupted - check data\n";
	throw e;
      }
      if (displevel>0) dispflag=true;
    }

    UMinLBFGS(const UMinLBFGS<Scalar> & umin)
      : ls(umin.ls), fx(umin.fx),
	firststep(umin.firststep),
	invhessianscale(umin.invhessianscale),
	maxupdates(umin.maxupdates),
	maxitn(umin.maxitn),
	tol(umin.tol),
	dispflag(umin.dispflag),
	str(umin.str) {}

    virtual ~UMinLBFGS() {}

    void setState( const Vector<Scalar> & s ) { 
      Vector<Scalar> & x = fx.getPoint();
      x.copy(s);
    }

    Vector<Scalar> & getState() {
      return fx.getPoint();
    }

    bool run() {
      try {

	CountingIterationTable<Scalar> 
	  stop(maxitn,
	       tol*fx.getGradient().norm(),
	       fx,
	       dispflag,
	       str);
	
	LBFGSStep<Scalar> step(ls,fx,
			       firststep,
			       invhessianscale,
			       maxupdates,
			       str);

	LoopAlg Umin(step, stop);
	return Umin.run();
      } 
      catch(RVLException & e) {
	e << "called from UMinLBFGS::run()\n";
	throw e;
      } 
      catch( std::exception & e) {
	RVLException es;
	es << "Exception caught in UMinLBFGS::run() with error message";
	es << e.what(); 
	throw e;
      }
      return true;
    }
  };
}

#endif
