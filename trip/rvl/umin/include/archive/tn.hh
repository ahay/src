// newtonalg.H
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

#ifndef __RVL_OPTALG
#define __RVL_OPTALG

#include "cgalg.hh"
#include "lnsrchBT.hh"
#include "conoptstep.hh"

namespace RVLUmin {
  using namespace RVLAlg;
  using namespace RVL;

  /** A general framework for Newton--like optimization methods  
      with globalization in terms of a variable step length.
      
      The state in this algorithm is the vector x
      which is the current best point.
  */

template<class Scalar>
class NewtonStep: public ConOptStep<Scalar> {
public:
  NewtonStep(Functional<Scalar> const & f, Vector<Scalar> const & x0)
    : ConOptStep<Scalar>(f, x0)
  {}

  ~NewtonStep() {}

  virtual bool run() {
    try {
      Vector<Scalar> dir(ConOptStep<Scalar>::getFunctionalEvaluation().getDomain(), true);
      bool cd = calcDir(dir);
      bool cs = calcStep(dir);
      return (cd && cs);      }
    } catch(RVLException & e) {
      e << "called from NewtonStep::run()\n";
      throw e;
    } catch( std::exception & e) {
      RVLException es;
      es << "Exception caught in NewtonStep::run() with error message";
      es << e.what(); 
      throw e;
    }
    return true;
  }

protected:
  // Methods which must/may be implemented by children
  //@{

  // Fill in the vector with the search direction
  virtual bool calcDir(Vector<Scalar> & dir) = 0;

  //@}

  // Methods which may be implemented by children
  //@{

  /** Find the step length in this search direction,
      and updates x.
      Default implementation always sets the step to 1.
  */
  virtual bool calcStep(Vector<Scalar> & dir) {
    Scalar step = 1.0;  
    (this->getTrialPoint()).copy(this->getBasePoint());
    (this->getTrialPoint()).linComb(step, dir);
    return true;
  }

  //@}  
};

/** This Truncated Newton Algorithm is described in
    Numerical Optimization by Nocedal & Wright.
    It is listed as Algorithm 6.1 on page 140.

    It uses a CGStep to find a search direction and
    then a LineSearchAlg to determine the step length.
*/
template <class Scalar>
class TruncatedNewtonStep: public NewtonStep<Scalar> {
public:
  /** Takes the functional to be minimized and
      a vector containing the current location x*/
  TruncatedNewtonStep(Functional<Scalar> & _f, const Vector<Scalar> & _x): 
    NewtonStep<Scalar>(_f,_x), bs(_f.getDomain())
  {}

  /** Takes a vector containing the starting location x0
      and a functional evaluation to be manipulated.
  */
  TruncatedNewtonStep(const Vector<Scalar> & _x0, FunctionalEvaluation<Scalar> & _fx): 
    NewtonStep<Scalar>(_x0, _fx), bs(_fx.getDomain())
  {}

  ~TruncatedNewtonStep() {}

protected:

  /** Find the step length in this search direction,
      and update x and fx.
      
      Returning a non-positive step length is 
      an algorithmic failure, and will cause the run()
      method to return false.

      Use a backtracking line search.
  */
  //  virtual bool calcStep(Vector<Scalar> & dir, Scalar & step) {
  virtual bool calcStep(Vector<Scalar> & dir) {
    bs.set(dir,ConOptStep<Scalar>::getFunctionalEvaluation());
    bool res = bs.run();
    //    step = bs.getStep();  
    return res;
  }

  /** Fill in the vector with the search direction
      Find a search direction using CG.
  */
  virtual bool calcDir(Vector<Scalar> & dir) {
    Scalar rho, tol;
    dir.zero();
    FunctionalEvaluation<Scalar> & fx = ConOptStep<Scalar>::getFunctionalEvaluation();
    Vector<Scalar> g(fx.getDomain());
    g.copy(fx.getGradient());
    g.negate();

    CGStep<Scalar> cgs(dir, fx.getHessian(), g, rho);

    /* We want tol = min(0.5, sqrt(g.norm())) * g.norm() */
    tol = g.norm();
    if( tol < 0.25 ) {
      tol = sqrt(tol)*tol;
    } else {
      tol = 0.5 * tol;
    }
    
    tol = tol *tol; // We will compare r.norm2() to tol^2

    MinTerminator<Scalar> t1(rho, tol);
    MinTerminator<Scalar> t2(cgs. getCurvature(), cgs.getCurvatureTol());
    OrTerminator t(t1,t2);

    LoopAlg la(cgs, t);

    return la.run();
  }
 
  BacktrackingLineSearchAlg<Scalar> bs;
};

}

#endif
