// gps.H
// created by ADP 12/22/03
// header file for the Generalized Pattern Search Algorithm

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

#ifndef __RVLALG_G_PATTERN_SEARCH_H 
#define __RVLALG_G_PATTERN_SEARCH_H  

#include "alg.hh"
#include "terminator.hh"
#include "CountElements.hh"
#include "identity.hh"
#include "functional.hh"
#include "conoptstep.hh"

#define inf numeric_limits<Scalar>::infinity()

namespace RVLUmin {
  using namespace RVLAlg;
  using RVL::Vector;
  using RVL::Space;
  using RVL::Functional;


  /** Abstract base class for a Search.  Adds step size access to ConOptStep */
template<class Scalar>
class LocalSearch : public ConOptStep<Scalar> {
private:
  bool success;
public:
  LocalSearch(Functional<Scalar> & _f, const Vector<Scalar> & startpoint)
    : ConOptStep<Scalar>(_f,startpoint), success(false) {}
  
  LocalSearch( FunctionalEvaluation<Scalar> & _fx, const Vector<Scalar> & startpoint)
    : ConOptStep<Scalar>(_fx, startpoint), success(false) {}

  virtual void setStepSize(Scalar step) = 0;
  virtual Scalar getStepSize() = 0;

  void setValue() { success=false; }
  void setValue(bool _success) { success=_success; }
  bool getValue() const { return success; }  
};

/** Abstract Base class for a local poll.  Implements run() and
    step size functions, but children are required to
    provide search directions.  As these search directions are not
    required to be fixed, or even have a constant number between runs,
    this should cover most Polls.
*/
template<class Scalar>
class LocalPoll : public LocalSearch<Scalar> {
protected:
  
  bool mini; // if true, minimize.  if false, maximize
  bool success;
  Scalar Delta;

  /** Return a reference to the ith search direction.  Since the 
      set of search directions at any given time is finite,
      this makes sense.
      
      Each search direction will be scaled when doing the linear combination,
      so the returned vectors should be on a reference scale with
      step size 1.  Thus, the norm of each should be roughly 1.

      valid inputs: \f$0 \leq i < n\f$ where n is the number returned
      by getNumSearchDir().
  */
  virtual Vector<Scalar> & getSearchDir(int i) = 0;

  /** Current number of search directions.  Must be a positive number  */
  virtual int getNumSearchDir() = 0;

public:
  LocalPoll(Functional<Scalar> & _f, Vector<Scalar> & startpoint, 
	    Scalar stepsize, bool minimize = true) 
    : LocalSearch<Scalar>(_f,startpoint), mini(minimize), Delta(stepsize), success(false);
  {}

  LocalPoll(FunctionalEvaluation<Scalar> & _fx, Vector<Scalar> & startpoint, 
	    Scalar stepsize, bool minimize = true) 
    : LocalSearch<Scalar>(_fx,startpoint), mini(minimize), Delta(stepsize), success(false);
  {}

  void setStepSize(Scalar step) {
    Delta = step;
  }
  virtual Scalar getStepSize() { return Delta; }

  void setValue() { success=false; }
  void setValue(bool _success) { success=_success; }
  bool getValue() const { return success; }

  /** Evaluate the function at each best+searchdir until an improvement is found
   */
  virtual run() {
    this->setValue();
    int i = 0;
    int n = getNumSearchDir();
    FunctionalEvaluation<Scalar> & fx = ConOptStep<Scalar>::getFunctionalEvaluation();
    Scalar fx0 = fx.getValue();
    Vector<Scalar> & x = this->getTrialPoint();
    const Vector<Scalar> & x0 = this->getBasePoint();  
    if( mini) {
      while(!(this->getValue()) && (i < n) ) {
	x.copy(x0);
	x.linComb(Delta, getSearchDir(i)); 
	if( fx.getValue() < fx0 ) {
	  this->setValue(true);
	}
	i++;
      }
    } else {
      while(!(this->getValue()) && (i < n) ) {
	x.copy(x0);
	x.linComb(Delta, getSearchDir(i)); 
	if( fx.getValue() > fx0 ) {
	  this->setValue(true);
	}
	i++;
      }
    }
  }

};

/**  This algorithm implements the two stage Generalized Pattern 
 Search step.  It first tries a global search.  If the search fails, it
 then performs a local poll.  If the poll fails, we have a mesh optimum
 and we return false.  A return value of true indicates that a better point
 was found.
*/
template <class Scalar>
class GPSStep : public LocalSearch<Scalar> {
protected:

  LocalSearch<Scalar> & searchstep;
  LocalPoll<Scalar> & pollstep;
  bool success;

public:

  /** NOTE:  Might be able to avoid passing F and startpoint.
      Seems silly to duplicate effort.
   */

  GPSStep( Functional<Scalar> & _f,
	   Vector<Scalar> & startpoint,
	   LocalSearch<Scalar> & search,
	   LocalPoll<Scalar>  & poll)
    : LocalSearch<Scalar>(_f, startpoint), searchstep(search), pollstep(poll) {}

  GPSStep( FunctionalEvaluation<Scalar> & _fx,
	   Vector<Scalar> & startpoint,
	   LocalSearch<Scalar> & search,
	   LocalPoll<Scalar>  & poll)
    : LocalSearch<Scalar>(_fx, startpoint), searchstep(search), pollstep(poll) {}

  void run() {
    searchstep.set(this->getBasePoint(),
		   this->getFunctionalEvaluation());
    searchstep.run();
    this->setValue(searchstep.getValue());
    if (!searchstep().getValue()) {
      // else poll      
      pollstep.set(this->getBasePoint(),
		   this->getFunctionalEvaluation());
      pollstep.run(); // return poll success
      this->setValue(pollstep.getValue());
    }
  }

  void setValue() { success=false; }
  void setValue(bool _success) { success=_success; }
  bool getValue() const { return success; }

  void setStepSize(Scalar step) {
    searchstep.setStepSize(step);
    pollstep.setStepSize(step);
  }
  virtual Scalar getStepSize() {
    return pollstep.getStepSize();
  }
};

  /** Bundles together GPSstep and some usual terminators.  A convienence class. */

template <class Scalar>
class GPSAlg : public LocalSearch<Scalar> {
protected:

  Scalar Delta;
  
  GPSStep<Scalar> step;
  MinTerminator<Scalar> mterm;
  SteppedIterationTable<Scalar> iterm;
  OrTerminator term;
  Scalar frac;

public:

  GPSAlg(Functional<Scalar> & _f, 
	 Vector<Scalar> & startpoint, 
	 LocalSearch<Scalar> & search,
	 LocalPoll<Scalar>  & poll,
	 Scalar tol = 1e-8,
	 Scalar fraction = 0.5,
	 Scalar initstepsize = 1
	 ) 
    : LocalSearch<Scalar>(_f,startpoint), Delta(initstepsize), step(_f,startpoint,search,poll), 
      mterm(Delta, tol), iterm(ConOptStep<Scalar>::getFunctionalEvaluation(), Delta), term(iterm, mterm),
      frac(fraction) 
  {}

  GPSAlg(FunctionalEvaluation<Scalar> & _fx, 
	 Vector<Scalar> & startpoint, 
	 LocalSearch<Scalar> & search,
	 LocalPoll<Scalar>  & poll,
	 Scalar tol = 1e-8,
	 Scalar fraction = 0.5,
	 Scalar initstepsize = 1
	 ) 
    : LocalSearch<Scalar>(_fx,startpoint), Delta(initstepsize), step(_fx,startpoint,search,poll), 
      mterm(Delta, tol), iterm(ConOptStep<Scalar>::getFunctionalEvaluation(), Delta), term(iterm, mterm),
      frac(fraction) 
  {}


  ~GPSAlg() {}

  void setStepSize(Scalar stepsize) {
    Delta = stepsize;
    step.setStepSize(stepsize);
  }
  virtual Scalar getStepSize() {
    return Delta;
  }

  void run() {
    Vector<Scalar> x0(this->getBasePoint());
    FunctionalEvaluation<Scalar> & fx = ConOptStep<Scalar>::getFunctionalEvaluation();
    Vector<Scalar> & x = this->getTrialPoint();
    while(!term.query()) {
      step.set(x0, fx);
      step.run();
      if( !step.getValue() ) {
	Delta = Delta*frac;
	step.setStepSize( Delta );
      } else {
	x0.copy(x);
      }
    }
  }
};

// Some concrete examples of search and poll steps

/** A simple compass search.  The ith search direction is the ith column of the
    identity matrix.
*/
template<class Scalar>
class CompassSearch : public LocalPoll<Scalar> {
protected:
  const Space<Scalar>& s;
  Vector<Scalar> lastsearchdir;
  int n;

  Vector<Scalar> & getSearchDir(int i) {
    if (i < n ) { 
      RVL::SetToColumnOfIdentity<Scalar> ei(i);
      lastsearchdir.eval(ei);
    } else {
      RVL::SetToColumnOfIdentity<Scalar> ei(i-n);
      lastsearchdir.eval(ei);
      lastsearchdir.negate();
    }
    return lastsearchdir;
  }

  int getNumSearchDir() {
    return n;
  }

  CompassSearch();
public:
  /** Takes a reference space in which the search will be done */
  CompassSearch( Functional<Scalar> & f, Vector<Scalar> & start, 
		 Scalar stepsize = 1, bool minimize = true)  
    : LocalPoll<Scalar>(f, start, stepsize, minimize), s(f.getDomain()), lastsearchdir(s) 
  {   
    CountElements<Scalar> c;
    lastsearchdir.eval(c);
    n = 2*c.getValue();
  }

};

  /** A placeholder search which doesn't do anything and never succeeds. */
template<class Scalar>
class DummySearch : public LocalSearch<Scalar> {
protected:
  Scalar Delta;

  DummySearch();

public:
  /** Search in the given space, starting from 0 */
  DummySearch( Functional<Scalar> & _f, Vector<Scalar> startpoint)
    : LocalSearch<Scalar>(_f,startpoint) {}

  void setStepSize(Scalar step) {
    Delta = step;
    
  }
  virtual Scalar getStepSize() {
    return Delta;
  }


  /** The DummySearch always fails */
  void run() {}

};

}
#endif
