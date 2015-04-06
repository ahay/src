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


#ifndef __ALG_CONOPTSTEP_H_
#define __ALG_CONOPTSTEP_H_

#include "alg.hh"
#include "functional.hh"

namespace RVLUmin {
  using namespace RVLAlg;
  using RVL::Functional;
  using RVL::Vector;

  /** Base class for Continuous Optimization Steps.
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
  class ConOptStep: public StateAlg<Vector<Scalar> > {
  private:
    Vector<Scalar> ownedx0; // the owned copy of the initial point, for use with setState
    Vector<Scalar> ownedx; // the owned copy of the trial point, for use with setState
    
    FunctionalEvaluation<Scalar> ownedfx; // the owned functional evaluation

    ConstDynRef< Vector<Scalar> > x0ref; // a reference to the initial point currently being used
    DynRef< Vector<Scalar> > xref; // a reference to the initial point currently being used
    DynRef< FunctionalEvaluation<Scalar> > fxref;// a reference to the initial point currently being used

    Scalar scale;
    Scalar gradscale;
    bool scaleset;

    // disable the copy constructor
    ConOptStep(const ConOptStep<Scalar> & cos);  

  public:
    /** Can be initialized with just a Functional.  Uses the
	domain of the functional to initialize workspace */

    ConOptStep( Functional<Scalar> & _f, 
		Scalar _scale = 0.0,
		Scalar _gradscale = 0.0
		)
      : ownedx0(_f.getDomain()), 
	ownedx(_f.getDomain()), ownedfx(_f,ownedx), 
	scale(0.0), gradscale(0.0), scaleset(false)
    {
      if( _scale != 0.0) {
	scale = _scale;
	scaleset = true;
      }
      if( _gradscale != 0.0) {
	gradscale = _gradscale;	
      }
    }

    /** Initialize with the functional and the starting point for the 
	optimization.
    */
    ConOptStep( Functional<Scalar> & _f, 
		const Vector<Scalar> & x0,
		Scalar _scale = 0.0,
		Scalar _gradscale = 0.0)
      : ownedx0(x0), ownedx(x0), ownedfx(_f,ownedx), scaleset(true) 
    {
	x0ref.set(ownedx0);
	xref.set(ownedx);
	fxref.set(ownedfx);
	
	if( _scale != 0.0) {
	  scale = _scale;
	} else {
	scale = ownedfx.getValue();
	}
	if( _gradscale != 0.0) {
	  gradscale = _gradscale;	
	} else {
	  gradscale = ownedfx.getGradient().norm();
	}
    }

    /** Constructor with the same effect as calling set() on these inputs.  Stores
	the necessary references.
    */
    ConOptStep( const Vector<Scalar> & x0, FunctionalEvaluation<Scalar> & fx,
		Scalar _scale = 0.0,
		Scalar _gradscale = 0.0)
      : ownedx(x0), ownedx0(x0), ownedfx(fx), scaleset(true)
    {
      x0ref.set(x0);
      xref.set(fx.getPoint());
      fxref.set(fx);

      if( _scale != 0.0) {
	scale = _scale;
      } else {
	scale = ownedfx.getValue();
      }
      if( _gradscale != 0.0) {
	gradscale = _gradscale;	
      } else {
	gradscale = ownedfx.getGradient().norm();
      }
    }

    virtual ~ConOptStep() {}

    /** Set up optimization problem with a point at which to start the 
	step, and a functional evaluation to modify as the optimization proceeds. 
    */
    void set(const Vector<Scalar> & x0,
	     FunctionalEvaluation<Scalar> & fx)
    {
      try {
	x0ref.set(x0);
	xref.set(fx.getPoint());
	fxref.set(fx);

	if(!scaleset) {
	  scale = fx.getValue();
	  gradscale = fx.getGradient().norm();
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from LineSearchAlg::set\n";
	throw e;
      }
    }

    /** Return a reference to the base point for the step. */
    const Vector<Scalar> & getBasePoint() { 
      try {
	return x0ref.get();
      }
      catch (RVLException & e) {
	e<<"\ncalled from LineSearchAlg::getBasePoint\n";
	throw e;
      }
    }

    /** Return a reference to the point at which the step terminated */
    Vector<Scalar> & getTrialPoint() {
      try {
	return xref.get();
      } catch (RVLException & e) {
	e<<"\ncalled from LineSearchAlg::getTrialPoint\n";
	throw e;
      }
    }

    /** Return a reference to the functional evaluation at the trial point */
    FunctionalEvaluation<Scalar> & getFunctionalEvaluation() {
      try {
	return fxref.get();
      }
      catch (RVLException & e) {
	e<<"\ncalled from ConOptStep::getFunctionalEvaluation\n";
	throw e;
      }
    }    
    
    Scalar getScale() {
      if(!scaleset) {
	scale = fxref.get().getValue();
	gradscale = fxref.get().getGradient().norm();
	scaleset = true;
      }
      return scale;
    }

    Scalar getGradScale() {
      if(!scaleset) {
	scale = fxref.get().getValue();
	gradscale = fxref.get().getGradient().norm();
	scaleset = true;
      }
      return gradscale;
    }

    /**  inherited from StateAlg.  Uses input as x0 and x.  sets references appropriately */
    virtual void setState(const Vector<Scalar> & x) {
      ownedx.copy(x);
      ownedx0.copy(x);
      xref.set(ownedx);
      x0ref.set(ownedx0);
      fxref.set(ownedfx);
    }

    /**  inherited from StateAlg.  Returns getTrialPoint() */
    virtual const Vector<Scalar> & getState() {
      return xref.get();
    }
  };

}


#endif
