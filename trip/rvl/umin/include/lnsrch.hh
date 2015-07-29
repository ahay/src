// lnsrch.H
// created by ADP
// last modified 06/17/04
// considerable modification 2007-2011 WWS

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


#ifndef __RVL_ALG_LINE_SEARCH_
#define __RVL_ALG_LINE_SEARCH_

#include "alg.hh"
#include "functional.hh"
#include "terminator.hh"
#include "table.hh"

namespace RVLUmin {

  using namespace RVL;
  using namespace RVLAlg;

  // forward declaration;
  template<typename Scalar>
  class LineSearchAlg;

  /**
     Base class for line search algorithms. Acts as a container for
     base point, gradient at base point, search direction, and
     functional evaluation at a current point. The run() method is
     deferred to child classes. The functional evaluation is a
     reference to an external object, which can be further manipulated
     by other program units - the purpose of the run() method is to
     change this object.

     On call to run,
     <ul>
     <li>x0 = base point for search</li>
     <li>g0 = gradient at base point</li>
     <li>dx = search direction</li>
     <li>fx = functional evaluation at x0</li>
     <li>step = initial (trial) step</li>
     </ul>

     On return from run,
     <ul>
     <li>x  = minimizer estimate</li>
     <li>fx = functional evaluation at x</li>
     </ul>

     getStep() returns the step taken.

     28.04.07: This all ought to be rewritten in Alg form. The LSBase
     should turn into an LSStep class, and LS should be rewritten to
     dynamically allocate a LSBase, and combine it with a terminator
     in a LoopAlg.
  */

  template <class Scalar>
  class LineSearchAlgBase: public Algorithm, public Terminator {

  private:

    Scalar f0;
    Vector<Scalar> x0;
    Vector<Scalar> g0;
    LineSearchAlg<Scalar> const & lsalg;

    LineSearchAlgBase();

  protected:

    FunctionalEvaluation<Scalar> & fx;
    Vector<Scalar> const & dx;

    // virtual copy constructor, in order that LSAlg can have real one
    virtual LineSearchAlgBase<Scalar> * clone() const = 0;

  public:

    LineSearchAlgBase(LineSearchAlg<Scalar> const & _lsalg,
		      FunctionalEvaluation<Scalar> & _fx,
		      Vector<Scalar> const & _dx)
	: 
	f0(_fx.getValue()),
	x0(_fx.getPoint()),
	g0(_fx.getGradient()),
	lsalg(_lsalg),
	fx(_fx), 
	dx(_dx) {} //cerr<<"from LineSearchAlgBase constructor\n";}
    LineSearchAlgBase(const LineSearchAlgBase<Scalar> & ls)
	: f0(ls.f0), x0(ls.x0), g0(ls.g0), lsalg(ls.lsalg), fx(ls.fx), dx(ls.dx) {}
    virtual ~LineSearchAlgBase() {}

    virtual ostream & write(ostream & str) const {
      str<<"lnsrch base class"<<endl;
      return str;
    }

    // note that Terminator::query and Algorithm::run are abstract

    virtual Scalar getStep() const = 0;
    Scalar getBaseValue() const { return f0; }
    Vector<Scalar> const & getBasePoint() const { 
      return x0; }
    Vector<Scalar> const & getBaseGradient() const { 
      return g0; }
    Vector<Scalar> & getTrialPoint() const { 
      return fx.getPoint(); }
    Vector<Scalar> const & getSearchDirection() const { 
      return dx; }
    FunctionalEvaluation<Scalar> & getFunctionalEvaluation() const { 
      return fx; }
    //    bool checkMinStep() const { return lsalg.checkMinStep(); }
    Scalar getMinStep() const { return lsalg.getMinStep(); }
  };

  /** Abstract handle class template for line searches. Instantiates
      its captive LineSearchAlgBase data member via protected pure
      virtual constructor (build(...)), and implements run() by
      delegation to this internal Alg. Also manages two key items:
      <ul> <li>Stores step last used, as scale for subsequent line
      searches - thus acts as persistent line search memory </li>
      <li>Minimum step length scale (relative to base point length),
      so that step can be checked to avoid destructive cancellation
      </li> </ul>

  */
      
  template<typename Scalar>
  class LineSearchAlg: public Algorithm, public Terminator {
  private:
    // the actual algorithm
    mutable LineSearchAlgBase<Scalar> * ls;
    // step used in last successful line search, initialized somehow
    mutable Scalar firststep;
    // minimum permitted step
    Scalar minsteptol;

  protected:
    /** Protected initialization - to be supplied by subclass. Note
	that LSAlgBase class is responsible for sensibly extracting an
	initial step from the firststep arg, which is initially the
	max Scalar, then the step actually taken in the last
	successful line search. */
    virtual LineSearchAlgBase<Scalar> * 
    build(FunctionalEvaluation<Scalar> & fx,
	  Vector<Scalar> const & dx,
	  Scalar firststep) = 0;

  public:
    /** constructor:
	@param _firststep = initial step length, modified by each run to be the last successful step
	@param _minsteptol = minimum permitted step
    */
    LineSearchAlg(Scalar _firststep = numeric_limits<Scalar>::max(),
		  Scalar _minsteptol = numeric_limits<Scalar>::epsilon())
      : ls(NULL), firststep(_firststep), minsteptol(_minsteptol) {}
    LineSearchAlg(LineSearchAlg<Scalar> const & lsh)
      : ls(NULL), firststep(lsh.firststep), minsteptol(lsh.minsteptol) {
      if (lsh.ls) ls=(lsh.ls)->clone();
    }
    ~LineSearchAlg() {
      if (ls) delete ls;
    }

    /** this method makes this class an abstract factory */
    void initialize(FunctionalEvaluation<Scalar> & fx,
		    Vector<Scalar> const & dx) {
      if (ls) delete ls;
      ls = build(fx,dx,firststep);
    }

    /** run whatever you've got, then record step actually taken */
    void run() {
      if (ls) { 
	ls->run();
	firststep=ls->getStep();
	return;
      }
      RVLException e;
      e<<"Error: LineSearchAlg::run\n";
      e<<"this alg not initialized\n";
      throw e;
    }

    /** query functions */
    bool query() { 
      if (ls) {
	return ls->query();
      }
      return false;
    }

    Scalar getStep() const { return ls->getStep(); }

    Scalar getBaseValue() const { 
      if (ls) return ls->getBaseValue();
      RVLException e;
      e<<"Error: LineSearchAlg::getBaseValue()\n";
      e<<"this alg not initialized\n";
      throw e;
    }

    Vector<Scalar> const & getBasePoint() const { 
      if (ls) return ls->getBasePoint(); 
      RVLException e;
      e<<"Error: LineSearchAlg::getBasePoint()\n";
      e<<"this alg not initialized\n";
      throw e;
    }

    Vector<Scalar> const & getBaseGradient() const { 
      if (ls) return ls->getBaseGradient(); 
      RVLException e;
      e<<"Error: LineSearchAlg::getBaseGradient()\n";
      e<<"this alg not initialized\n";
      throw e;
    }

    Vector<Scalar> & getTrialPoint() const { 
      if (ls) return ls->getTrialPoint(); 
      RVLException e;
      e<<"Error: LineSearchAlg::getTrialPoint()\n";
      e<<"this alg not initialized\n";
      throw e;
    }

    Vector<Scalar> const & getSearchDirection() const { 
      if (ls) return ls->getSearchDirection(); 
      RVLException e;
      e<<"Error: LineSearchAlg::getSearchDirection()\n";
      e<<"this alg not initialized\n";
      throw e;
    }

    FunctionalEvaluation<Scalar> & getFunctionalEvaluation() const {
      if (ls) return ls->getFunctionalEvaluation(); 
      RVLException e;
      e<<"Error: LineSearchAlg::getFunctionalEvaluation()\n";
      e<<"this alg not initialized\n";
      throw e;
    }

    /** Checks that the step is not senselessly small. Size of step
	times direction vector is measured either relative to the
	length of the current iterate, if the latter is nonzero, or
	absolutely if the current iterate happens to be the zero
	vector. In either case returns true if stepsize is greater
        than this measure, else false. */
    /*
    bool checkMinStep() const {
      try {
	Scalar foo = this->getStep();
	return (foo > minsteptol);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LineSearchAlg::checkMinStep\n";
	throw e;
      }
    }
    */
    Scalar getMinStep() const { return minsteptol; }

  };

}

#endif
