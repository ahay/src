// conoptstep.H
// created by ADP 5/27/04

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

#ifndef __ALG_UMINSTEP_H_
#define __ALG_UMINSTEP_H_

#include "lnsrch.hh"
#include "alg.hh"
#include "functional.hh"

namespace RVLUmin {

  using namespace RVLAlg;
  using RVL::Functional;
  using RVL::Vector;
  using RVL::Vector;
  using RVL::FunctionalEvaluation;

  /** Abstract interface for computation of search directions, in
      support of descent methods.  Since computation can succeed or
      fail, must have character of terminator.
  */

  template<typename Scalar>
  class UMinDir: public Terminator {
  public:
    UMinDir() {}
    UMinDir(UMinDir<Scalar> const &) {}
    virtual ~UMinDir() {}

    // note that Terminator::query must be defined by instantiable subclass 

    /** Returns search direction in mutable first argument */
    virtual void calcDir(Vector<Scalar> & dir,
			 FunctionalEvaluation<Scalar> & fx) = 0;

    /** Use data generated during line search to update any internals */
    virtual void updateDir(LineSearchAlg<Scalar> const & ls) = 0;

    /** Reinitialize direction computation */
    virtual void resetDir() = 0;

    /** verbose output stream */
    virtual ostream & write(ostream & str) const = 0;
  };

  /** Base class for Unconstrained Minimization step algorithms
      with globalization via line search.
.
      Provides the added functionality of maintaining the
      functional evaluations and states correctly when
      set, and providing external access to such states.
      The functional evaluation may be set to avoid recomputing
      result in the calling object.

      On call,
      x0 = base point for search
      fx = functional evaluation at x0
        
      On return,
      x  = minimizer estimate
      fx = functional evaluation at x

      The state is always the point x.
  */
  template<class Scalar>
  class UMinStepLS: public Algorithm, public Terminator {
  private:
    
    FunctionalEvaluation<Scalar> & fx; 
    UMinDir<Scalar> & dc;
    LineSearchAlg<Scalar> & ls;
    bool ans;
    ostream & str;

    /** Find the step length in this search direction,
	and updates x.
	Default implementation always sets the step to 1.
	
	Returning a non-positive step length is 
	an algorithmic failure, and will cause the run()
	method to return false.
    */

    virtual void calcStep(Vector<Scalar> & dir) {
      try {
	//Line Search
	bool tried_steepest_descent = false;
	// set up and run line search
	ls.initialize(fx,dir);
      	ls.run();
	bool res = ls.query();
	if (!res) { 
	  dc.updateDir(ls);
	  ans = ans || dc.query(); 
	}
	// line search failed, try steepest descent
	else {
	  if (!tried_steepest_descent) {
	    str<<"UMinStep: attempting steepest descent restart\n";
	    // This resets the evaluation point
	    Vector<Scalar> & x = fx.getPoint();
	    x.copy(ls.getBasePoint());
	    dir.copy(ls.getBaseGradient());
	    dir.negate();
	    dc.resetDir();
	    // this costs a redundant gradient computation!!
	    ls.initialize(fx,dir);
	    ls.run();
	    if (!ls.query()) { 
	      dc.updateDir(ls);
	      ans= ans || dc.query();
	    }
	    else {
	      ans = true;
	    }
	    tried_steepest_descent=true;
	  }
	  else {
	    ans = true;
	  }
	}

      } catch(RVLException & e) {
	e << "\ncalled from UMinStepLS::calcStep()\n";
	throw e;
      }
    }

  public:

    /** Initialize with the functional and the starting point for the 
	optimization.
    */
    UMinStepLS(FunctionalEvaluation<Scalar> & _fx, 
	       UMinDir<Scalar> & _dc,
	       LineSearchAlg<Scalar> & _ls,
	       ostream & _str = cout)
      : fx(_fx), dc(_dc), ls(_ls), ans(false), str(_str) {}

    UMinStepLS(const UMinStepLS<Scalar> & cos)
      : fx(cos.fx), dc(cos.dc), ls(cos.ls), ans(cos.ans), str(cos.str) {}

    virtual ~UMinStepLS() {}

    /** Return a reference to the base point for the step. */
    Vector<Scalar> const & getBasePoint() { return ls.getBasePoint(); }

    /** Return a reference to the gradient at the bast point. */
    Vector<Scalar> const & getBaseGradient() { return ls.getBaseGradient(); }

    /** Return a reference to the functional evaluation at the trial point */
    FunctionalEvaluation<Scalar> & getFunctionalEvaluation() { return fx; }
    
    void run() {
      try {
	Vector<Scalar> dir(fx.getDomain(), true);
	dc.calcDir(dir,fx);
	calcStep(dir);
      } 
      catch(RVLException & e) {
	e << "called from UMinStepLS::run()\n";
	throw e;
      } 
      catch( std::exception & e) {
	RVLException es;
	es << "Exception caught in UMinStepLS::run() with error message";
	es << e.what(); 
	throw e;
      }
    }

    bool query() { return ans; }

  };

}


#endif
