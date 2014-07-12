#ifndef __RVLALG_LINFITSM_L2_H
#define __RVLALG_LINFITSM_L2_H

/** Given an Operator F and a Vector d in the range of op,
    implements the function
    \f$$
    f(x) = \inf_{dx} \|DF(x)dx - d\|^2
    \f$$
    as an RVL::Functional. The linear least squares solver is specified by 
    policy.
*/

#include "alg.hh"
#include "terminator.hh"
#include "linop.hh"
#include "table.hh"

using namespace RVLAlg;

namespace RVLUmin {

  using namespace RVL;
  using namespace RVLAlg;    

  template<typename Scalar, typename LSPolicy, typename LSPolicyData> 
  class LinFitLSSM: public Functional<Scalar>, public LSPolicy {
    
    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  private:
    
    Operator<Scalar> const & op;        // operator  
    LinearOp<Scalar> const & preop;     // preconditioner
    LinearOp<Scalar> const & helmop;    // smoothing op applied to gradient
    Vector<Scalar> const & d;           // data 
    Vector<Scalar> const & x0;          // input initial linear solution

    bool refine;                        // refine as in Kern & Symes 1994
    mutable  Vector<Scalar> dx;         // preimage of linear solution
    mutable  Vector<Scalar> dltx;       // linear solution
    mutable bool applied;
    ostream & str;
    
  protected:

    void apply(const Vector<Scalar> & x, 
	       Scalar & val) const {
      try {
	/*         if (applied) {
		   RVLException e;
		   e<<"Error: LinFitLS::apply(x,val)\n";
		   e<<"already applied, may not alter\n";
		   throw e;
		   }
	*/
	atype rnorm;
	atype nrnorm;
	// access Operator through OperatorEvaluation
	OperatorEvaluation<Scalar> opeval(op,x);

	// Get Derivative of Operator
	LinearOp<Scalar> const & lop = opeval.getDeriv();

	// Composition of lop and preop
	OpComp<Scalar> gop(preop,lop);

	Vector<Scalar> tmp(gop.getDomain());
	tmp.zero();

        // for given initial solution
    Vector<Scalar> d0(lop.getRange());
    lop.applyOp(x0,d0);
    d0.linComb(1.0,d,-1.0);

	dx.zero();

	// build least square solver , solve for dx
	OperatorEvaluation<Scalar> gopeval(gop,x0);

	Algorithm * solver 
	  = LSPolicy::build(dx,gopeval.getDeriv(),d0,rnorm,nrnorm,str);

	solver->run();

	// get the value of objective function
	val = 0.5*rnorm*rnorm;

	preop.applyOp(dx,dltx);
        dltx.linComb(1.0,x0);
      
	applied = true;
	delete solver;
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinFitLSSM::apply\n";
	throw e;
      }
    } 

    void applyGradient(const Vector<Scalar> & x, 
		       Vector<Scalar> & g) const {
        try{
	  if(!applied){
            Scalar val;
            this->apply(x,val);
	  }
	  OperatorEvaluation<Scalar> opeval(op,x);
	  LinearOp<Scalar> const & lop = opeval.getDeriv();
	  SymmetricBilinearOp<Scalar> const & sblop = opeval.getDeriv2();
	  Vector<Scalar> dltd(lop.getRange());
	  Vector<Scalar> gtmp(g.getSpace());
	  // compute dltx and dltd = DF * dltx - d
	  lop.applyOp(dltx,dltd);
	  dltd.linComb(-1.0,d);
	  // naive computation of gradient
	  sblop.applyAdjOp(dltx,dltd,gtmp);
	  
	  // compute and add correction term to gradient
	  if (refine) {
	    atype rnorm;
	    atype nrnorm;
	    OpComp<Scalar> gop(preop,lop);
	    Vector<Scalar> tmp(gop.getDomain());	 
	    Vector<Scalar> dx1(gop.getDomain());
	    tmp.zero();
	    dx1.zero();
	    OperatorEvaluation<Scalar> gopeval(gop,tmp);
	    // solve DF * dx = dltd in LS sense 
	    Algorithm * solver = LSPolicy::build(dx1,gopeval.getDeriv(),dltd,rnorm,nrnorm,str);
	    solver->run();
	    delete solver;
	  
	    Vector<Scalar> tmp2(g.getSpace());
	    Vector<Scalar> dx2(preop.getRange());
	    preop.applyOp(dx1,dx2);
	    // compute and add correction term tmp to gradient g
	    sblop.applyAdjOp(dx2,d,tmp2);
	    gtmp.linComb(1.0, tmp2);
	  }
          cerr << "\n before helmop! \n" << endl;
          helmop.applyOp(gtmp,g);
          cerr << "\n after helmop! \n" << endl;
        }
        catch (RVLException & e) {
	  e<<"\ncalled from LinFitLSSM::applyGradient\n";
            throw e;
        }
        
    }

    void applyHessian(const Vector<Scalar> & x,
		      const Vector<Scalar> & dx, 
		      Vector<Scalar> & dy) const {}

    Functional<Scalar> * clone() const {
      return new LinFitLSSM<Scalar,LSPolicy,LSPolicyData>(*this);
    }

  public:

    /* typical policy data 
	       atype _rtol,
	       atype _nrtol,
	       int _maxcount,
    */
    LinFitLSSM(Operator<Scalar> const & _op,
	     LinearOp<Scalar> const & _preop,
	     LinearOp<Scalar> const & _helmop,
	     Vector<Scalar> const & _d,
             Vector<Scalar> const & _x0,
	     LSPolicyData const & s,
	     bool _refine=false,
	     ostream & _str=cerr)
	: LSPolicy(), op(_op), helmop(_helmop), preop(_preop), d(_d), x0(_x0), 	  
	  refine(_refine), dx(preop.getDomain()), dltx(preop.getRange()), 
	  applied(false), str(_str) {
      try{
	dx.zero();
	LSPolicy::assign(s);
	if (s.verbose) {
	  str<<"\n";
	  str<<"==============================================\n";
	  str<<"LinFitLSSM constructor - ls policy data = \n";
	  s.write(str);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinFitLSSM::Constructor\n";
	throw e;
      }
    }

    LinFitLSSM(LinFitLSSM<Scalar,LSPolicy,LSPolicyData> const & f) 
	: LSPolicy(f), op(f.op), preop(f.preop), helmop(f.helmop), d(f.d), x0(f.x0),refine(f.refine),
	  dx(f.dx), dltx(f.dltx), applied(f.applied), str(f.str) {}

    const Space<Scalar> & getDomain() const { return op.getDomain(); }

    Scalar getMaxStep(const Vector<Scalar> & x,
		      const Vector<Scalar> & dx) const {
      try {
	return op.getMaxStep(x,dx);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinFitLSSM::getMaxStep\n";
	throw e;
      }
    }

    Vector<Scalar> const & getLSSoln() const { return dltx; }

    ostream & write(ostream & str) const {
      str<<"LinFitLSSM: \n";
      str<<"*** operator:\n";
      op.write(str);
      str<<"*** data vector:\n";
      d.write(str);
      return str;
    }
  };
}
#endif
