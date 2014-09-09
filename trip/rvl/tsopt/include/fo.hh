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
#ifndef __RVL_STATEFO
#define __RVL_STATEFO

#include "alg.hh"
#include "dtstep.hh"
#include "local.hh"
#include "contentpackage.hh"

namespace TSOpt {

  using RVL::RVLException;
  using RVL::LocalDataContainer;
  using RVL::UnaryLocalFunctionObject;
  using RVL::UnaryLocalFunctionObjectScalarRedn;

  /** Copies \ref RVL::LocalDataContainer to state and control arrays of RnState. */
  template<typename State>
  class LDCtoStateFOR: public UnaryLocalFunctionObjectScalarRedn<float,bool> {
  private:
    State & s;
    LDCtoStateFOR();
    LDCtoStateFOR(LDCtoStateFOR<State> const &);
  public:
    LDCtoStateFOR(State & _s)
      : s(_s) {}
    ~LDCtoStateFOR() {}

    void operator() (LocalDataContainer<float> const & x) {
      try {
	// initialize dimension - it=0
	// <ME?> what if s already exists?
	if (!( s.getData() && s.getSize()>=x.getSize())) {
	  RVLException e;
	  e<<"Error: TSOpt::LDCtoStateFOR::operator()\n";
	  e<<"state data not initialized, or too small to accomodate input\n";
	  e<<"state size="<<s.getSize()<<" input size="<<x.getSize()<<"\n";
	  throw e;
	}
	memcpy(s.getData(),x.getData(),x.getSize()*sizeof(float));
      }
      catch (RVLException & e) {
	e<<"\ncalled from LDCtoStateFOR::operator()\n";
	throw e;
      } 
    }

    string getName() const { return "LDCtoStateFOR"; }
  };

  /** Copies \ref RVL::LocalDataContainer from state array of RnState. */
  template<typename State>
  class StatetoLDCFO: public UnaryLocalFunctionObject<float> {
  private:
    State const & s;
    StatetoLDCFO();
    StatetoLDCFO(StatetoLDCFO<State> const &);
  public:
    StatetoLDCFO(State const & _s)
      : s(_s) {}
    ~StatetoLDCFO() {}

    void operator() (LocalDataContainer<float> & x) {
      try {
	// check dimension
	if (s.getSize() != x.getSize()) {
	  RVLException e;
	  e<<"Error: StatetoLDCFO::operator()\n";
	  e<<"external, interal output dims don't agree\n";
 
	  throw e;
	}
	memcpy(x.getData(),s.getData(),x.getSize()*sizeof(float));
      }
      catch (RVLException & e) {
	e<<"\ncalled from StatetoLDCFO::operator()\n";
	throw e;
      } 
    }

    string getName() const { return "StatetoLDCFO"; }
  };

}

#endif
