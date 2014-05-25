#ifndef __RVLALG_LINFIT_L2_H
#define __RVLALG_LINFIT_L2_H

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
  class LinFitLS: public Functional<Scalar>, public LSPolicy {
    
    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  private:
    
    Operator<Scalar> const & op;        // operator  
    LinearOp<Scalar> const & preop;     // preconditioner
    Vector<Scalar> const & d;           // data 
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

	Vector<Scalar> x0(gop.getDomain());

	x0.zero();

	dx.zero();

	// build least square solver , solve for dx
	OperatorEvaluation<Scalar> gopeval(gop,x0);

	Algorithm * solver 
	  = LSPolicy::build(dx,gopeval.getDeriv(),d,rnorm,nrnorm,str);

	solver->run();

	// get the value of objective function
	val = 0.5*rnorm*rnorm;

	preop.applyOp(dx,dltx);
      
	applied = true;
	delete solver;
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinFitLS::apply\n";
	throw e;
      }
    } 

    void applyGradient(const Vector<Scalar> & x, 
		       Vector<Scalar> & g) const {
        try{
	  if(!applied){
            Scalar val;
            this->apply(x,val);
            cerr << "\n val=" << val << endl;
	  }
	  OperatorEvaluation<Scalar> opeval(op,x);
	  LinearOp<Scalar> const & lop = opeval.getDeriv();
	  SymmetricBilinearOp<Scalar> const & sblop = opeval.getDeriv2();
	  Vector<Scalar> dltd(lop.getRange());
	  // compute dltx and dltd = DF * dltx - d
	  lop.applyOp(dltx,dltd);
          //  cerr << "\n dltx.norm()=" << dltx.norm() << endl;
	  dltd.linComb(-1.0,d);
          //cerr << "\n dltd.norm()=" << dltd.norm() << endl;
	  // naive computation of gradient
	  sblop.applyAdjOp(dltx,dltd,g);
          //  cerr << "\n g.norm()=" << g.norm() << endl;
	  
	  // compute and add correction term to gradient
          cerr << "\n LinFitLS refine=" << refine << endl;
	  if (refine) {
          cerr << "\n LinFitLS: inside refine branch\n" ;
	    atype rnorm;
	    atype nrnorm;
	    OpComp<Scalar> gop(preop,lop);
	    Vector<Scalar> x0(gop.getDomain());	 
	    Vector<Scalar> dx1(gop.getDomain());
	    x0.zero();
	    dx1.zero();
	    OperatorEvaluation<Scalar> gopeval(gop,x0);
	    // solve DF * dx = dltd in LS sense 
	    Algorithm * solver = LSPolicy::build(dx1,gopeval.getDeriv(),dltd,rnorm,nrnorm,str);
	    solver->run();
	    delete solver;
	  
	    Vector<Scalar> tmp(g.getSpace());
	    Vector<Scalar> dx2(preop.getRange());
	    preop.applyOp(dx1,dx2);
	    // compute and add correction term tmp to gradient g
	    sblop.applyAdjOp(dx2,d,tmp);
            cerr << "\n LinFitLS:applyGradient  term1.norm = " << g.norm() << endl;
            cerr << "\n LinFitLS:applyGradient  term2.norm = " << tmp.norm() << endl;
	    g.linComb(1.0, tmp);
            cerr << "\n LinFitLS:applyGradient  g.norm = " << g.norm() << endl;
	  }
        }
        catch (RVLException & e) {
	  e<<"\ncalled from LinFitLS::applyGradient\n";
            throw e;
        }
        
    }

    void applyHessian(const Vector<Scalar> & x,
		      const Vector<Scalar> & dx, 
		      Vector<Scalar> & dy) const {}

    Functional<Scalar> * clone() const {
      return new LinFitLS<Scalar,LSPolicy,LSPolicyData>(*this);
    }

  public:

    /* typical policy data 
	       atype _rtol,
	       atype _nrtol,
	       int _maxcount,
    */
    LinFitLS(Operator<Scalar> const & _op,
	     LinearOp<Scalar> const & _preop,
	     Vector<Scalar> const & _d,
	     LSPolicyData const & s,
	     bool _refine=false,
	     ostream & _str=cerr)
	: LSPolicy(), op(_op), preop(_preop), d(_d), 	  
	  refine(_refine), dx(preop.getDomain()), dltx(preop.getRange()), 
	  applied(false), str(_str) {
      try{
	dx.zero();
	LSPolicy::assign(s);
	if (s.verbose) {
	  str<<"\n";
	  str<<"==============================================\n";
	  str<<"LinFitLS constructor - ls policy data = \n";
	  s.write(str);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinFitLS::Constructor\n";
	throw e;
      }
    }

    LinFitLS(LinFitLS<Scalar,LSPolicy,LSPolicyData> const & f) 
	: LSPolicy(f), op(f.op), preop(f.preop), d(f.d), refine(f.refine),
	  dx(f.dx), dltx(f.dltx), applied(f.applied), str(f.str) {}

    const Space<Scalar> & getDomain() const { return op.getDomain(); }

    Scalar getMaxStep(const Vector<Scalar> & x,
		      const Vector<Scalar> & dx) const {
      try {
	return op.getMaxStep(x,dx);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LinFitLS::getMaxStep\n";
	throw e;
      }
    }

    Vector<Scalar> const & getLSSoln() const { return dltx; }

    ostream & write(ostream & str) const {
      str<<"LinFitLS: \n";
      str<<"*** operator:\n";
      op.write(str);
      str<<"*** data vector:\n";
      d.write(str);
      return str;
    }
  };
}
#endif
