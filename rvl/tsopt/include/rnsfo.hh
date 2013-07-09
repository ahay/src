/*************************************************************************

Copyright Rice University, 2004, 2005, 2006, 2007, 2008
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
#ifndef __RVL_RNSTATEFO
#define __RVL_RNSTATEFO

#include "rn.hh"
#include "rns.hh"
#include "realrns.hh"
#include "alg.hh"
#include "dtstep.hh"
#include "contentpackage.hh"
#include "local.hh"

namespace TSOpt {

  using RVL::RVLException;
  using RVL::LocalDataContainer;
  using RVL::UnaryLocalFunctionObject;
  using RVL::UnaryLocalFunctionObjectConstEval;

  /** Copies \ref RVL::LocalDataContainer to state and control arrays of RnState. */
  class LDCtoRnFOR: public UnaryLocalFunctionObjectConstEval<float> {
  private:
    RnState & s;
    LDCtoRnFOR();
    LDCtoRnFOR(LDCtoRnFOR const &);
  public:
    LDCtoRnFOR(RnState & _s)
      : s(_s) {}
    ~LDCtoRnFOR() {}

    using RVL::UnaryLocalConstEval<float>::operator();
    void operator() (LocalDataContainer<float> const & x);

    string getName() const { return "LDCtoRnFOR"; }
  };

  /** Copies \ref RVL::LocalDataContainer from state array of RnState. */
  class RntoLDCFO: public UnaryLocalFunctionObject<float> {
  private:
    RnState const & s;
    mutable float * out;
    RntoLDCFO();
    RntoLDCFO(RntoLDCFO const &);
  public:
    RntoLDCFO(RnState const & _s)
      : s(_s), out(NULL) {}
    ~RntoLDCFO();

    using RVL::UnaryLocalEvaluation<float>::operator();
    void operator() (LocalDataContainer<float> & x);
    string getName() const { return "RntoLDCFO"; }
  };

  /** Copies \ref RVL::LocalDataContainer to state and control arrays of RealRnState. */
  class LDCtoRealRnFOR: public UnaryLocalFunctionObjectConstEval<double> {
  private:
    RealRnState & s; 
    double * c;
    LDCtoRealRnFOR();
    LDCtoRealRnFOR(LDCtoRealRnFOR const &);
  public:
    LDCtoRealRnFOR(RealRnState & _s)
      : s(_s), c(NULL) 
    { 
      c = new double(s.getrealrn().nc);
      memcpy(c,s.getrealrn().c, s.getrealrn().nc*sizeof(double));  
    }

    ~LDCtoRealRnFOR() {delete c;}

    using RVL::UnaryLocalConstEval<double>::operator();    
    void operator() (LocalDataContainer<double> const & x) {

      try {
	// initialize dimension - it=0
	// <ME?> what if s already exists?
	if (!( s.getrealrn().u && s.getrealrn().c && s.getrealrn().nu>=1 && s.getrealrn().nc>=1))  {
	  s.initialize(x.getSize(),x.getSize());
	}

	size_t nc = s.getrealrn().nc;
	
	memcpy(s.getrealrn().c, c, nc*sizeof(double));
	memcpy(s.getrealrn().u, x.getData(), x.getSize()*sizeof(double));
	
	RealTime t;
	// <ME?> Is this incorrect? Should be at the initial time of 
	//       the simulation, whether fwd or adj
	t = 0.0;
	s.setTime(t);
	
	cout << "x.getData()[0] = " << x.getData()[0] << endl;	
	cout << "s= " << endl;
	s.write(cout);

	
      }
      catch (RVLException & e) {
	e<<"\ncalled from LDCtoRealRnFOR::operator()\n";
	throw e;
      } 
      

    }

    string getName() const { return "LDCtoRealRnFOR"; }
  };

  /** Copies \ref RVL::LocalDataContainer from state array of RealRnState. */
  class RealRntoLDCFO: public UnaryLocalFunctionObject<double> {
  private:
    RealRnState const & s;
    mutable double * out;
    RealRntoLDCFO();
    RealRntoLDCFO(RealRntoLDCFO const &);
  public:
    RealRntoLDCFO(RealRnState const & _s)
      : s(_s), out(NULL) { }
    ~RealRntoLDCFO();

    using RVL::UnaryLocalEvaluation<double>::operator();
    void operator() (LocalDataContainer<double> & x);
    string getName() const { return "RealRntoLDCFO"; }
  };


}

#endif
