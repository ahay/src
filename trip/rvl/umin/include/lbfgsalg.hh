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

#include "uminstep.hh"
#include "linop.hh"


namespace RVLUmin{

  using namespace RVLAlg;
  using RVL::Vector;
  using RVL::LinearOp;
  using RVL::Functional;
  using RVL::FunctionalEvaluation;
  using RVL::Table;
  
  /** LMBFGSOp implements the limited memory BFGS approximation
      to the inverse Hessian of a twice-differentiable function.  This
      approximation uses local changes to the gradient to gradually
      build up an estimate of the Hessian for use in nonlinear
      optimization problems.
    
      For details of the algorithm, see the paper
    
      "Updating Quasi-Newton Matrices with Limited Storage" by Jorge Nocedal,
      Math. of Computation, Vol. 35, no. 151, p.p. 773--782.
    
      Note that the operator is an approximation to the <i>inverse</i> Hessian,
      so the the apply method computes the quasi-Newton step.
    
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

  protected:
  
    LinearOp<Scalar> * clone() const {
      return new LBFGSOp<Scalar>(*this);
    }

    // apply computes the image of the operator on x, giving y.
    void apply(const Vector<Scalar> & x,
	       Vector<Scalar> & y) const {
      try {
      
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
	this->applyOp(y,x);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LBFGSOp::applyAdj\n";
	throw e;
      }
    }

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
    ~LBFGSOp() {}

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
    void update(const Vector<Scalar> & x,
		const Vector<Scalar> & xnext,
		const Vector<Scalar> & g,
		const Vector<Scalar> & gnext) {
      try {
	if (MaxNum==0) return;
      
	if( CurAllocated < MaxNum ) {
	  CurAllocated++;
	}

	S[CurNum].copy(xnext);
	S[CurNum].linComb(-1.0,x);
	Y[CurNum].copy(gnext);
	Y[CurNum].linComb(-1.0,g);
	
	if (ProtectedDivision<Scalar>
	    (1.0,S[CurNum].inner(Y[CurNum]),rho[CurNum])) {
	  RVLException e;
	  e<<"LBFGSOp::update\n";
	  e<<"zerodivide in first protected div\n";
	  e<<"either secant vector vanishes, or change in gradient vector vanishes, or \n";
	  e<<"secant is perpindicular to change in gradient\n";
	  e<<"----- row/col index = "<<CurNum<<"\n";
	  e<<"----- secant vector:\n";
	  S[CurNum].write(e);
	  e<<"----- delta grad: \n";
	  Y[CurNum].write(e);
	  //	  throw e;
	  e<<"revert to initial inv Hessian approximation\n";
	  reset();
	}

	Scalar tmp;
	if (ProtectedDivision<Scalar>
	    (1.0,rho[CurNum]*(Y[CurNum].normsq()),tmp)) {
	  RVLException e;
	  e<<"LBFGSOp::update\n";
	  e<<"zerodivide in second protected div\n";
	  e<<"change in gradient vanishes\n";
	  //	  throw e;
	  e<<"revert to initial inv Hessian approximation\n";
	  reset();
	}

	//	(ProtectedDivision<Scalar>
	//	 (1.0,rho[CurNum]*(Y[CurNum].normsq()),tmp));

	H0.setScale(tmp);
	CurNum++;
	if(CurNum == MaxNum) CurNum = 0;

      }
      catch (RVLException & e) {
	e<<"\ncalled from LBFGSOp::update\n";
	throw e;
      }
    }

    // reset sets the operator to the initial inverse Hessian approximation.
    void reset() { CurNum = 0; CurAllocated = 0; }

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
  class LBFGSDir : public UMinDir<Scalar> {

  private:

    LBFGSDir();

    LBFGSOp<Scalar> H;
    bool ans;
    ostream & str;

  public:

    // Fill in the vector with the search direction
    // nothing should go wrong here, if H has been updated.
    void calcDir(Vector<Scalar> & dir,
		 FunctionalEvaluation<Scalar> & fx) {
      try {
	// why isn't H defined to be -H?
	const Vector<Scalar> & grad = fx.getGradient();
	H.applyOp(grad,dir);
	dir.negate();
      } catch(RVLException & e) {
	e << "called from LBFGSStep::calcDir()\n";
	throw e;
      }
    }

    void updateDir(LineSearchAlg<Scalar> const & ls) {
      try {
	H.update(ls.getBasePoint(),
		 ls.getTrialPoint(),
		 ls.getBaseGradient(),
		 ls.getFunctionalEvaluation().getGradient());
      }
      catch (RVLException & e) {
	e<<"\ncalled from LBFGSDir::updateDir\n";
	throw e;
      }
    }

    void resetDir() { H.reset(); }
    
    LBFGSDir(Space<Scalar> const & dom,
	     Scalar InvHessianScale = 1.0,
	     int MaxUpdates = 5,
	     ostream & _str=cout)
      :  UMinDir<Scalar>(),	
	 H(dom, 
	   InvHessianScale, 
	   MaxUpdates),  ans(false), str(_str) {}
  
    LBFGSDir(LBFGSDir<Scalar> const & x) 
      : UMinDir<Scalar>(x), H(x.H) {}
    ~LBFGSDir() {} 

    bool query() { return ans; }

    ostream & write(ostream & str) const {
      str<<"LBFGSDir\n";
      //  H.write(str);
      return str;
    }
  };

}

#endif
