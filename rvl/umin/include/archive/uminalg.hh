
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

#ifndef __ALG_UMINALG_H_
#define __ALG_UMINALG_H_

#include "uminstep.hh"

namespace RVLUmin {

  using namespace RVLAlg;
  using RVL::Functional;
  using RVL::Vector;

  template<class Scalar> 
  class UMinAlg : public StateAlg<Vector<Scalar> > {
  public:
    virtual void run() = 0;

    virtual Vector<Scalar> & getState()  = 0;
    virtual Vector<Scalar> const & getState() const  = 0;

  };

  /* temp disabled
  
  template<class Scalar>
  class UMinAlgTR : public UMinAlg<Scalar> {
  private:

    Functional<Scalar> const & f; // objective function 
    UMinDir<Scalar> & dc;         // direction update method
    Vector<Scalar> & x;           // initial estimate on call, solution on return

    // Terminator parameters
    int maxitn;
    Scalar tol;
    int displevel;
    bool dispflag;
    ostream & str;
    
    UMinAlgTR();

  public:
    UMinAlgTR(Functional<Scalar> const & _f,
	      Vector<Scalar> & _x,
	      UMinDir<Scalar> & _dc,
	      int _maxitn = 20,
	      Scalar _tol = 0.01,
	      bool _dispflag = false,
	      ostream & _str = cout
	      )
      : f(_f), x(_x), dc(_dc),
	maxitn(_maxitn),
	tol(_tol),
	dispflag(_dispflag),
        str(_str) {}

    UMinAlgTR(UMinAlgTR<Scalar> const & umin)
      : f(umin.f), x(umin.x), dc(umin.dc),
	maxitn(umin.maxitn),
	tol(umin.tol),
	dispflag(umin.dispflag),
	str(umin.str) {}
	       
    virtual ~UMinAlgTR() {}

    void run() {
      try {
	
	FunctionalEvaluation<Scalar> fx(f,x);
      
	CountingIterationTable<Scalar> 
	  stop(maxitn,
	       //	       tol*fx.getGradient().norm(),
	       tol,
	       fx,
	       dispflag,
	       str);
	
	UMinStepTR<Scalar> step(fx,dc,str);
	LoopAlg Umin(step, stop);
	Umin.run();
      } 
      catch(RVLException & e) {
	e << "called from UMinAlgTR::run()\n";
	throw e;
      } 
      catch( std::exception & e) {
	RVLException es;
	es << "Exception caught in UMinAlgTR::run() with error message";
	es << e.what(); 
	throw e;
      }
      return true;
    } 

    Vector<Scalar> & getState() { return x; }
    Vector<Scalar> const & getState() const { return x; }

    ostream & write(ostream & str) const {
      str<<"UMinAlgTR: unconstrained minimization with TR-CG algorithm\n";
      str<<"Objective function:\n";
      f.write(str);
      str<<"Solution vector:\n";
      x.write(str);
      return str;
    }


  };



  template<class Scalar>
  class UMinAlgLS : public UMinAlg<Scalar> {
  private:

    Functional<Scalar> const & f; // objective function 
    Vector<Scalar> & x;           // initial estimate on call, solution on return
    UMinDir<Scalar> & dc;         // direction update method
    LineSearchAlg<Scalar> & ls;   // line search method

    // Terminator parameters
    int maxitn;
    Scalar tol;
    int displevel;
    bool dispflag;
    ostream & str;
    
    UMinAlgLS();

  public:

    UMinAlgLS(Functional<Scalar> const & _f,
	    Vector<Scalar> & _x,
	    UMinDir<Scalar> & _dc,
	    LineSearchAlg<Scalar> & _ls,
	    int _maxitn = 20,
	    Scalar _tol = 0.01,
	    bool _dispflag = false,
	    ostream & _str = cout
	    )
      : f(_f), x(_x), dc(_dc), ls(_ls),
	maxitn(_maxitn),
	tol(_tol),
	dispflag(_dispflag),
        str(_str) {}

    UMinAlgLS(UMinAlgLS<Scalar> const & umin)
      : f(umin.f), x(umin.x), dc(umin.dc), ls(umin.ls),
	maxitn(umin.maxitn),
	tol(umin.tol),
	dispflag(umin.dispflag),
	str(umin.str) {}
	       
    virtual ~UMinAlgLS() {}

    Vector<Scalar> & getState() { return x; }
    Vector<Scalar> const & getState() const { return x; }

    void run() {
      try {

	FunctionalEvaluation<Scalar> fx(f,x);
	
	CountingIterationTable<Scalar> 
	  stop(maxitn,
	       //	       tol*fx.getGradient().norm(),
	       tol,
	       fx,
	       dispflag,
	       str);

	
	UMinStepLS<Scalar> step(fx,dc,ls,str);

	LoopAlg Umin(step, stop);

	Umin.run();
      } 
      catch(RVLException & e) {
	e << "called from UMinAlgLS::run()\n";
	throw e;
      } 
      catch( std::exception & e) {
	RVLException es;
	es << "Exception caught in UMinAlgLS::run() with error message";
	es << e.what(); 
	throw e;
      }
      return true;
    }

    ostream & write(ostream & str) const {
      str<<"UMinAlgLS: unconstrained minimization with line search globalization\n";
      str<<"Objective function:\n";
      f.write(str);
      str<<"Solution vector:\n";
      x.write(str);
      str<<"Direction update:\n";
      dc.write(str);
      str<<"Line Search algorithm:\n";
      ls.write(str);
      return str;
    }
  };
}

#endif
