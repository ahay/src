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

#ifndef __RVL_ALG_NEWTON_H_
#define __RVL_ALG_NEWTON_H_

#include "alg.hh"
#include "op.hh"
#include "linop.hh"
#include "linearsolver.hh"

namespace RVLAlg {

  using RVL::Vector;
  using RVL::OperatorEvaluation;

  /** Classic Newton solver step for finding x so that \f$F(x)=0\f$.  Requires 
      that the derivative of F be invertable and F know how to evaluate
      such an inverse.
  */
template<class Scalar>
class NewtonSolverStep: public StateAlg<Vector<Scalar> > {
protected:
  Vector<Scalar> x;
  OperatorEvaluation<Scalar> eval;
  OperatorEvaluation<Scalar> & opeval;
  Vector<Scalar> & p;

public:
  /** Build a solver step starting at the zero vector in the domain */
  NewtonSolverStep( RVL::OperatorWithInvertibleDeriv<Scalar> & op ) 
    : x(op.getDomain(), true), eval(op, x), opeval(eval), p(x){}

  /** Build a solver step using the given evaluation, allowing for non-zero
      initial guesses
  */
  NewtonSolverStep( OperatorEvaluation<Scalar> & _opeval ) 
    : x(_opeval.getPoint()), eval(_opeval), 
      opeval(_opeval), p(opeval.getPoint()) {
    // Need to check for invertable derivative?    
  }

  ~NewtonSolverStep() {}

  void run() { 
    try {

      Vector<Scalar> dx(opeval.getDomain());

      try {
	const Vector<Scalar> & Fx = opeval.getValue();
	const RVL::LinearOpWithInverse<Scalar> & DF = 
	  dynamic_cast<const RVL::LinearOpWithInverse<Scalar> &>
	  (opeval.getDeriv());
	DF.getInvOp().apply(Fx, dx);
      } catch (bad_cast) {
	RVL::RVLException e;
	e << "Bad cast in NewtonSolverStep::run() - Derivative is not invertable\n";
	throw e;
      }
      p.linComb(-1.0, dx);
      return;
    } 
    catch( RVL::RVLException & e ) {
      e << "in NewtonSolverStep::run()";
      throw e;
    }
  }

  /*  virtual void setState(const Vector<Scalar> & st) { x.copy(st); }*/
  virtual Vector<Scalar> & getState() { return x;}
  virtual const Vector<Scalar> & getState() const { return x;}
};

  /** A more generic Newton--Krylov method for finding x so that\f$F(x)=0\f$.
      The LinearSolver is used to find solutions of $DF(x) * dx = F(x)$.
      This allows the solving of systems for which the derivative is
      invertable, but the operator does not necessarily know how to find
      such an inverse.
  */

template<class Scalar>
class NewtonKrylovSolverStep: public StateAlg<Vector<Scalar> >, public Terminator {
protected:
  Vector<Scalar> x;
  OperatorEvaluation<Scalar> eval;
  OperatorEvaluation<Scalar> & opeval;
  LinearSolver<Scalar> & linsolve;
  Vector<Scalar> & p;
  
  mutable bool ans;

public:
  /** Build a solver step starting at the zero vector in the domain */
  NewtonKrylovSolverStep( RVL::Operator<Scalar> & op, LinearSolver<Scalar> & _linsolve ) 
    : x(op.getDomain(), true), eval(op, x), opeval(eval), 
      linsolve(_linsolve), p(opeval.getPoint()), ans(false) 
  {}

  /** Build a solver step using the given evaluation, allowing for non-zero
      initial guesses
  */
  NewtonKrylovSolverStep( OperatorEvaluation<Scalar> & _opeval, LinearSolver<Scalar> & _linsolve ) 
    : x(_opeval.getPoint()), eval(_opeval), opeval(_opeval), 
      linsolve(_linsolve), p(opeval.getPoint()), ans(false)
  {}

  ~NewtonKrylovSolverStep() {}

  bool query() { return ans; }

  void run() { 
    try {
      Vector<Scalar> dx(opeval.getDomain());
      linsolve.setSystem(opeval.getDeriv(), dx, opeval.getValue());
      linsolve.run();
      if (linsolve(query)) {
	p.linComb(-1.0, dx);      
	ans=true;
	return;
      }
      else {
	ans=false;
	return;
      }
    } 
    catch( RVL::RVLException & e ) {
      e << "in NewtonSolverStep::run()";
      throw e;
    }
  }

  /*  virtual void setState(const Vector<Scalar> & st) { p.copy(st); }*/
  virtual Vector<Scalar> & getState() { return p; }  
  virtual const Vector<Scalar> const & getState() { return p; }

};

}

#endif
