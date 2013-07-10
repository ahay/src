#ifndef __RVL_ALG_LINEAR_SOLVER_H_
#define __RVL_ALG_LINEAR_SOLVER_H_

#include "linop.hh"
#include "alg.hh"
#include "cgalg.hh"
#include "ioterm.hh"

namespace RVLAlg {

  using namespace RVLUmin;
  using RVL::Vector;
  using RVL::LinearOp;
 
  /** A simple interface for linear solvers.  Allows the system to be
      set dynamically. In effect, this is a factory class, which does
      not expose its product.

      An iterative solver of a linear system does not generally define
      a linear, nor even differentiable, map. For this reason, this
      type is not a subtype of any of the vector function classes, but
      is merely an algorithm.

      The terminator inheritance allows subclasses to internalize the
      success/failure evaluation. Apps will use it to determine
      whether the iteration converged, or the system was solved to
      adequate precision.
  */
  template<class Scalar>
  class LinearSolver : public Algorithm, public Terminator {

  public:     

    virtual ~LinearSolver() {}

    // really a build method 
    virtual void setSystem(LinearOp<Scalar> const & op,
			   Vector<Scalar> & sol, 
			   Vector<Scalar> const & rhs) = 0;

    virtual bool query() = 0;
    virtual void run() = 0;
  };

  /** Linear solver for symmetric systems, implemented via the CG
      algorithm */

  template<class Scalar>
  class CGLinearSolver : public LinearSolver<Scalar> {

    typedef typename ScalarFieldTraits<Scalar>::AbsType atype;

  private:
    SimpleCGAlg<Scalar> * myCG; 
    int maxIter;
    atype tol;
    atype maxstep;               
    ostream & str;
    mutable bool ans;

  public:

    CGLinearSolver(int _maxIter=10, atype _tol=1.e-8, atype _maxstep = numeric_limits<atype>::max(), ostream & _str = cout)	   
      : maxIter(_maxIter), tol(_tol), maxstep(_maxstep), str(_str), ans(false), myCG(NULL) {}
 
    CGLinearSolver(CGLinearSolver<Scalar> const & sv)
      : maxIter(sv.maxIter), tol(sv.tol), str(sv.str), ans(sv.ans), myCG(NULL) {}

    ~CGLinearSolver() { if (myCG) delete myCG; }

    void setSystem(LinearOp<Scalar> const & op,
		   Vector<Scalar> & sol, 
		   Vector<Scalar> const & rhs) {
      if (myCG) delete myCG;
      myCG = new SimpleCGAlg<Scalar>(sol,op,rhs,tol,maxIter,maxstep,str); 
    }

    bool query() { return ans; }
    void run() {
      try {
	if (!myCG) {
	  RVLException e;
	  e<<"Error: CGLinearSolver::run()\n";
	  e<<"linear system not initialized\n";
	  throw e;
	}
	myCG->run();
	ans=(myCG->getResidNormSq() < tol);
      }
      catch (RVLException & e) {
	e<<"\ncalled from CGLinearSolver::run()\n";
	throw e;
      }
    }
    
  };

}

#endif
