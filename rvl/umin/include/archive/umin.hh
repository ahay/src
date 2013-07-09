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

#ifndef __RVLALG_UMIN_CHOOSER_H
#define __RVLALG_UMIN_CHOOSER_H

#include "umintable.hh"
#include "lbfgsalg.hh"
#include "lnsrchBT.hh"
#include "uminstep.hh"
//#include "TRalg.hh" under construction

/** "Chooser" class for UMin methods - closed handle to StateAlg */

namespace RVLUmin {
  using namespace RVLAlg;
  using RVL::Vector;
  using RVL::Functional;
  using RVL::FunctionalEvaluation;

  template<class Scalar>
  class UMinMethod: public StateAlg< Vector<Scalar> > {
    
  private:

    Functional<Scalar> const & f;
    Vector<Scalar> & x;
    UMinTable<Scalar> & t;
    mutable int count;
    ostream & str;

    UMinMethod();
    UMinMethod(const UMinMethod<Scalar> &);

  public:

    UMinMethod(Functional<Scalar> const & _f,
	       Vector<Scalar> & _x,
	       UMinTable<Scalar> & _t,
	       ostream & _str = cout)
      : f(_f), x(_x), t(_t), str(_str), count(0) {}

    void run() {
      try {

	string upd;
	// check that update type is specified
	if (t.getValue("UpdateAlg",upd)) {
	  RVLException e;
	  e<<"Error: UMinMethod constructor\n";
	  e<<"update method not specified (key UpdatedAlg)\n";
	  throw e;
	}
	if (upd == "LBFGS") {
	  // BFGS update rule - also a terminator
	  LBFGSDir<Scalar> dir(f.getDomain(),
			       getValFromTable<Scalar>(t,"InvHessianScale"),
			       getValFromTable<int>(t,"MaxUpdates"),
			       str); 
	  // for LBFGS must specify line search method
	  // if (t.getValue("LineSearchAlg",lns)) {
	  //	    RVLException e;
	  //	    e<<"Error: UMinMethod constructor\n";
	  //	    e<<"line search method not specified (key LineSearchAlg) for LBFGS update\n";
	  //	    throw e;
	  //	  }
	  //if (lns == "DS") {
	  //    meth = new UMinLBFGS_DS<Scalar>(f,x,t,str);
	  //  }
	  //else
	  //	  if (lns == "BT") {

	  BacktrackingLineSearchAlg<Scalar>
	    ls(getValFromTable<int>(t,"MaxSample"),
	       getValFromTable<bool>(t,"DispFlag"),
	       str,
	       getValFromTable<Scalar>(t,"FirstStep"),
	       getValFromTable<Scalar>(t,"FcnDecreaseTol"),
	       getValFromTable<Scalar>(t,"StepReductionFactor"),
	       getValFromTable<Scalar>(t,"FractionOfMaxStep")
	       );

	  FunctionalEvaluation<Scalar> fx(f,x);

	  CountingIterationTable<Scalar> 
	    stop1(fx,
		 getValFromTable<int>(t,"MaxItn"),
		 getValFromTable<Scalar>(t,"AbsGradTol"),
		 getValFromTable<Scalar>(t,"RelGradTol"),
		 getValFromTable<bool>(t,"DispFlag"),
		 str);
	  
	  OrTerminator stop(dir,stop1);

	  UMinStepLS<Scalar> step(fx,dir,ls,str);

	  LoopAlg Umin(step, stop);

	  Umin.run();
	  
	  count=stop1.getCount();

	}
	else if (upd == "TRCG") {
	  RVLException e;
	  e<<"Error UMinMethod constructor\n";
	  e<<"TRCG currently not available\n";
	  throw e;
	  /*
	  cout << "Update = TRCG!" << endl;
	 
	  dir = new TRCGDir<Scalar>(f.getDomain(), NULL, 
				    getValFromTable<Scalar>(t,"InvHessianScale"), 
				    getValFromTable<int>(t,"MaxUpdates")); 
	 
	  meth = new UMinAlgTR<Scalar>
	    (f,x, *dir,
	     getValFromTable<int>(t,"MaxItn"),
	     getValFromTable<Scalar>(t,"GradTol"),
	     getValFromTable<bool>(t,"DispFlag"),
	     str);
	  */
	}
	else {
	  RVLException e;
	  e<<"Error: UMinMethod constructor\n";
	  e<<"UpdateAlg = "<<upd<<" has illegal value\n";
	  e<<"see docs for UMinMethod for currently available choices\n";
	  throw e;
	}	  
      }
      catch (RVLException & e) {
	e<<"\ncalled from UMinMethod constructor\n";
	throw e;
      }
    }

    ~UMinMethod() {}

    Vector<Scalar> & getState() {
      return x;
    }
    Vector<Scalar> const & getState() const {
      return x;
    }

    int getCount() const { return count; }
  };

}
#endif
