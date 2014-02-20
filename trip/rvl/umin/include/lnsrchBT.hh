//lnsrhcBT.H
// created by WWS
// last modified 06/09/04

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

#ifndef __RVL_ALG_LINE_SEARCH_BT
#define __RVL_ALG_LINE_SEARCH_BT

#include "table.hh"
#include "lnsrch.hh"

namespace RVLUmin {
  using namespace RVLAlg;
  using namespace RVL;

  /** Does a backtracking line search starting from a prescribed step,
      passed as argument firststep to the constructor. The initial
      step must be prescribed externally, as no scale information is
      available at the outset. Starts with backtrack loop: at each
      unsuccessful iteration, step <- step*gamma, with gamma <
      1. Stops when sufficient decrease condition is satisfied. If
      this happens on first step, goes into internal doubling loop: at
      each doubling iteration, step <- step*stretch, with 1 < stretch
      < gamma (last inequality to avoid cycling). Internal doubling
      stops if function value increases.

  */

  template <class Scalar>
  class BacktrackingLineSearchAlgBase: public LineSearchAlgBase<Scalar> {

  private:

    //    Scalar gamma; replaced by gamma1
    Scalar gamma1;
    //    Scalar con; replaced by eta1
    Scalar eta1;
    Scalar eta2;
    bool DispFlag;
    Scalar fudge;
    //    Scalar stretch; replaced by gamma2
    Scalar gamma2;
    int maxsteps;
    mutable Scalar step; // may be stretched on return
    bool ans;

    BacktrackingLineSearchAlgBase();

    ostream & str;

  protected:

    LineSearchAlgBase<Scalar> * clone() const {
      return new BacktrackingLineSearchAlgBase(*this);
    }

  public:

    /**
      Constructor arguments:
      @param lsalg const reference to line search algorithm
      @param fx current FunctionalEvaluation on call, updated on return
      @param dx const reference to search direction
      @param _step initial step on call, last successful step on return
      @param _eta1 lower G-A parameter > 0
      @param _eta2 upper G-A parameter > eta1
      @param _gamma1 trust region reduction factor < 1
      @param _gamma2 trust region expansion factor > 1, gamma1*gamma2 < 1 
      @param _DispFlag verbosity flag
      @param _fudge fraction of max step to permit
      @param _maxsteps max number of line search steps
      @param _str verbose output unit
    */
    BacktrackingLineSearchAlgBase
    (LineSearchAlg<Scalar> const & lsalg,
     FunctionalEvaluation<Scalar> & fx,
     Vector<Scalar> const & dx,
     Scalar _step,
     Scalar _eta1,
     Scalar _eta2,
     Scalar _gamma1,
     Scalar _gamma2,
     bool _DispFlag,
     Scalar _fudge,
     int _maxsteps,
     ostream & _str
     )
      : LineSearchAlgBase<Scalar>(lsalg,fx,dx),
	step(_step),
	gamma1(_gamma1),
	gamma2(_gamma2),
	eta1(_eta1),
	eta2(_eta2),
	DispFlag(_DispFlag),
	fudge(_fudge),
	maxsteps(_maxsteps),
	str(_str),
	ans(false) {}

    BacktrackingLineSearchAlgBase
    (const BacktrackingLineSearchAlgBase<Scalar> & ls)
      : LineSearchAlgBase<Scalar>(ls),
	step(ls.step),
        gamma1(ls.gamma1),
        gamma2(ls.gamma2),
        eta1(ls.eta1),
        eta2(ls.eta2),
	fudge(ls.fudge),
	maxsteps(ls.maxsteps),
        DispFlag(ls.DispFlag),
	str(ls.str), ans(ls.ans) {}

    virtual ~BacktrackingLineSearchAlgBase() {}
    
    ostream & write(ostream & str) const {
      str<<"bt\n";
      return str;
    }

    bool query() { return ans; }

    Scalar getStep() const { return step; }
    
    /**
       Checks for sufficient decrease, and while not found, shrink step
    */

    void run() {

      try {

	const Vector<Scalar> & x0 = this->getBasePoint();
	Vector<Scalar> const & dx = this->getSearchDirection();
	FunctionalEvaluation<Scalar> & fx = this->LineSearchAlgBase<Scalar>::getFunctionalEvaluation();
	Vector<Scalar> & x = fx.getPoint();
	Scalar minstep = this->getMinStep();

	//	cerr<<"fval in lnsrchBT::run"<<endl;

	Scalar fval = fx.getValue();
	Scalar gfdx = dx.inner(fx.getGradient());
	Scalar dxnorm = dx.norm();
	Scalar rate = gfdx/dxnorm;
	Scalar maxstp = fx.getMaxStep(dx);

	Scalar fruit = fudge*maxstp;
	step = (step < fruit)?step:fruit;

	if (DispFlag) {
	  str<<"BacktrackingLineSearchAlg::run:\n";
	  str<<"  initial function value = "<<fval<<"\n";
	  str<<"  initial descent rate   = "<<rate<<"\n";
	  str<<"  estimated step         = "<<step<<"\n";
	  str<<"  step vector norm       = "<<dxnorm<<"\n";
	  str<<"  comparison G-A value   = "<<fval+eta1*step*gfdx<<"\n";
	  str<<"  max feasible step      = "<<maxstp<<"\n";
	  str.flush();
	}

 	if (rate > 0.0) {
	  if (DispFlag) {
	    str<<"BacktrackingLineSearchAlg::run:\n";
	    str<<"  SEARCH DIRECTION IS NOT A DESCENT DIRECTION\n";
	    str.flush();
	  }
	  ans=true;
	  return;
	}

	//	cout<<"fruit = "<<fruit<<" step = "<<step<<endl;

	// check that step is not too small rel length of direction vector
        if (step<minstep*dxnorm) {
	  if (DispFlag) {
	    str<<"BacktrackingLineSearchAlg::run:\n";
	    str<<"  proposed step is too small rel search vector length\n";
	    str<<"  line search aborted\n";
	    str.flush();
	  }
	  ans=true;
	  return;
	}

	x.copy(x0);
	x.linComb(step, dx);
	
	int bt = 0;
	if (DispFlag) {
	  bt++;
	  str<<"BacktrackingLineSearchAlg::run:\n";
	  str<<"  backtrack iter "<<bt<<"\n";
	  str<<"  trial function value = "<<fx.getValue()<<"\n";
	  str<<"  trial step           = "<<step<<"\n";
	  str<<"  Min G-A overestimate = "<<fval+eta1*step*gfdx<<"\n";
	  str.flush();
	}        

	//	cout<<"LS: first step\n";
	//	x.write(cout);

	// while not sufficient decrease, shrink step 
	while( (fx.getValue() > fval + eta1*step*gfdx) && 
	       (step*dxnorm > minstep) && 
	       bt <= maxsteps) {
	  // if the new value is bigger than the old (esp. if it's much bigger),
	  // then a more nuanced estimate is sensible - use quadratic interpolation
	  if (fx.getValue() > fval && bt<2) {
	    Scalar tmpstep = -(gfdx*step*step)/(fx.getValue()-fval-step*gfdx);
	    str<<"  trial quadr. bt step = "<<tmpstep<<"\n";
	    if (tmpstep < minstep*dxnorm) {
	      step = step * gamma1;
	      if (DispFlag) {
		str<<"  quadratic bt step    = "<<tmpstep<<" too small, replace with\n";
		str<<"  linear bt step       = "<<step<<"\n";
		str.flush();
	      }
	    }
	    else {
	      step = tmpstep;
	      if (DispFlag) {
		str<<"  quadratic bt step    = "<<step<<" accepted\n";
		str.flush();
	      }
	    }
	  }
	  else {
	    step = step * gamma1;
	  }
	  x.copy(x0);
	  if (step*dxnorm < minstep) {
	    if (DispFlag) {
	      str<<"BacktrackingLineSearchAlg::run:\n";
	      str<<"  proposed step length = "<<step*dxnorm<<" is smaller than\n";
	      str<<"  minimum permitted = "<<minstep<<"\n";
	      str<<"  --- line search aborted\n";
	      str.flush();
	    }
	    ans=true;
	    return;
	  }	    
	  x.linComb(step, dx);

	  if (DispFlag) {
	    bt++;
	    str<<"BacktrackingLineSearchAlg::run:\n";
	    str<<"  backtrack iter "<<bt<<"\n";
	    str<<"  trial function value = "<<fx.getValue()<<"\n";
	    str<<"  Min G-A overestimate = "<<fval+eta1*step*gfdx<<"\n";
	    str<<"  trial step           = "<<step<<"\n";
	    str.flush();
	  }
	}

	// if we get to here, we've backtracked the max number of steps without
	// satisfying sufficient decrease
	// return "true", i.e. stop line search alg
	if (fx.getValue() > fval + eta1*step*gfdx && 
	    (step*dxnorm > minstep)) {
	  if (DispFlag) {
	    str<<"BacktrackingLineSearchAlg::run:\n";
	    str<<"  G-A criterion not satisfied, step not accepted\n";
	    str.flush();
	  }
	  if ( bt > maxsteps+1 ) { 
	    if (DispFlag) {
	      str<<"  Termination: maximum step count exceeded\n";
	      str.flush();
	    }
	  }
	  x.copy(x0);
	  ans=true;
	}
	// otherwise, line search has succeeded, check if we can do better
	else {
	  // internal doubling: lengthen step if success on first try - 
	  // requires additional workspace:
	  if( bt<2 ) {
	    if (DispFlag) {
	      str<<"BacktrackingLineSearchAlg::run:\n";
	      str<<"  Successful step\n";
	      str.flush();
	    }
	    Vector<Scalar> y(x.getSpace());  // vector in which to save last trial point
	    y.copy(x); 
	    Scalar stepsave = step;          // last step
	    Scalar fvalnext=fx.getValue();   // last function value
	    bt=0;                            // count reset for internal doubling
	    while( (fx.getValue() <= fval+eta2*step*gfdx) &&
		   (fx.getValue() <= fvalnext) &&
		   bt <= maxsteps) {
	      str<<"  successful internal doubling - try another\n";
	      stepsave=step;
	      y.copy(x);
	      step=step*gamma2;
	      fvalnext=fx.getValue();
	      x.copy(x0);
	      x.linComb(step, dx);
	      bt++;
	      if (DispFlag) {
		str<<"  internal doubling step "<<bt<<"\n";
		str<<"  trial function value = "<<fx.getValue()<<"\n";
		str<<"  Max G-A overestimate = "<<fval+eta2*step*gfdx<<"\n";	  
		str<<"  trial step           = "<<step<<"\n";
		str.flush();
	      }
	    }
	    // Stop doubling. If no improvement, backtrack.
	    if  ((fx.getValue() > fvalnext) && bt > 0) { 
	      step=stepsave;
	      x.copy(y);
	      if (DispFlag) {
		str<<"  last internal doubling step "<<bt<<" failed - backtrack\n";
	      }
	    }
	    else {
	      if (DispFlag) {
		str<<"  internal doubling terminated at step "<<bt<<"\n";
	      }		
	    }
	    if (DispFlag) {
	      str<<"  final function value = "<<fx.getValue()<<"\n";
	      str<<"  final step           = "<<step<<"\n";
	      str<<"  *** end line search ***\n";
	      str.flush();
	    }
	  }
	  ans=false;
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from BacktrackingLineSearchAlg::run\n";
	throw e;
      }
    }

  };

  /** Factory class for BacktrackingLineSearchAlgBase implementation
      of backtracking line search. Subclasses RVLUmin::LineSearchAlg
      abstract factory: build method accepts RVL::FunctionalEvaluation
      (function evaluated at current estimate of solution),
      RVL::Vector (step), both mutable references, and first step
      Scalar, and returns pointer to dynamically allocated
      BacktrackingLineSearchAlgBase object.

      Data members (initialized in constructor) define a backtracking line search:
      <ul>
      <li>firststep = initial step, stored by base class (RVLUmin::LineSearchAlg) and updated to last successful step with each run of line search</li> 
      <li>minsteptol = minimum permitted step, stored by base class (RVLUmin::LineSearchAlg) and used to sanity-test step. Typical value = 0.01</li>
      <li>eta1 = lower G-A parameter: if actual reduction in objective is less than eta1 * (linear prediction of reduction), then no update, step *= gamma1 ("backtracking"); else, step is (provisionally) accepted, pending internal doubling test (next bullet). Typical value = 0.01</li>
      <li>eta2 = upper G-A parameter: if actual reduction in objective is greater than eta2 * (linear prediction of reduction) AND this occurs on first line search step, then (a) step *= gamma2, (b) update recomputed, until either max line search iterates reached or eta2 test fails. At end of this "internal doubling" loop, step is accepted if eta1 test passes, else previous internal doubling update is accepted. Typical value = 0.9</li>  
      <li>gamma1 = step decrease factor. Typical value = 0.5</li>
      <li>gamma2 = step increase factor. Typical value = 1.8</li>
      <li>DispFlag = verbosity flag. Typical value = true</li>
      <li>fudge = max fraction of max step to permit. Typical value = 0.9</li>
      <li>maxsteps = max number of line search steps. Typical value = 10</li>
      <li>str = verbose output stream</li>
      </ul>

      Constraints on parameters:
      <ul>
      <li> 0 < eta1 < eta2 < 1</li>
      <li> 0 < gamma1 < 1 < gamma2</li>
      <li> gamma1*gamma2<1 (to prevent cycling)</li>
      <li> 0 < firststep < numeric_limits<Scalar>::max()</li>
      </ul>

      Typical values: eta1=0.01, eta2=0.5, gamma1=0.5, gamma2=1.8, firststep=1.0, fudge=0.9, maxsteps=10.

      The problematic choice is firststep - it depends on scales in
      the problem which are opaque from the point of view of this
      algorithm. In principle, the algorithm can recover from a bad
      choice of scale by either backtracking or internal doubling, but
      in practice this process can consume an indordinate number of
      function evaluations. The initial value of firststep should
      reflect whatever is known about the probable distance between
      the initial estimate of the solution and the actual optimizer.

   */
  template<typename Scalar>
  class BacktrackingLineSearchAlg: public LineSearchAlg<Scalar> {

  private:

    Scalar eta1;
    Scalar eta2;
    Scalar gamma1;
    Scalar gamma2;
    Scalar fudge;
    bool DispFlag;
    int maxsteps;
    ostream & str;

  protected:

    virtual LineSearchAlgBase<Scalar> * 
    build(FunctionalEvaluation<Scalar> & fx,
	  Vector<Scalar> const & dx,
	  Scalar firststep) {

      return new BacktrackingLineSearchAlgBase<Scalar>(*this,
						       fx,dx,firststep,
						       eta1,eta2,gamma1,gamma2,DispFlag,
						       fudge,maxsteps,str);
    }
    
  public:

    /**
      Constructor arguments:
      @param firststep initial step, stored by RVLUmin::LineSearchAlg base class
      @param minsteptol minimum permitted step, stored by RVLUmin::LineSearchAlg base class
      @param _eta1 lower G-A parameter > 0
      @param _eta2 upper G-A parameter > eta1
      @param _gamma1 trust region reduction factor < 1
      @param _gamma2 trust region expansion factor > 1, gamma1*gamma2 < 1 
      @param _DispFlag verbosity flag
      @param _fudge fraction of max step to permit
      @param _maxsteps max number of line search steps
      @param _str verbose output unit
    */
    BacktrackingLineSearchAlg(int _maxsteps=10,
			      bool _DispFlag = false,
			      Scalar firststep=1.0,
			      Scalar minsteptol=numeric_limits<Scalar>::min(),
			      Scalar _eta1=0.01,
			      Scalar _eta2=0.5,
			      Scalar _gamma1=0.5,
			      Scalar _gamma2=1.8,
			      Scalar _fudge=0.9,
			      ostream & _str = cout)
      : LineSearchAlg<Scalar>(firststep,minsteptol),
	DispFlag(_DispFlag),
	eta1(_eta1),
	eta2(_eta2),
	gamma1(_gamma1),
	gamma2(_gamma2),
	fudge(_fudge),
	maxsteps(_maxsteps),
        str(_str) {}

    /** copy constructor */
    BacktrackingLineSearchAlg(BacktrackingLineSearchAlg<Scalar> const & bt)
      : LineSearchAlg<Scalar>(bt),
	DispFlag(bt.DispFlag),
	eta1(bt.eta1),
	eta2(bt.eta2),
	gamma1(bt.gamma1),
	gamma2(bt.gamma2),
	fudge(bt.fudge),
	maxsteps(bt.maxsteps),
        str(bt.str) {}

    ~BacktrackingLineSearchAlg() {}

  };
}

#endif
