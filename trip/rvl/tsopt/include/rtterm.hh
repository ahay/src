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
#ifndef __RVL_RTTERM
#define __RVL_RTTERM

#include "rt.hh"
#include "tterm.hh"

namespace TSOpt {

  using RVL::RVLException;
  using RVLAlg::Terminator;

  /** FwdTimeTerm fully implemented with RealTime data member */
  template<typename State>
  class FwdRTimeTerm: public FwdTimeTerm<State> {
  private:
    RealTime t0;
    FwdRTimeTerm();
  public:
    FwdRTimeTerm(State const & t)
      : FwdTimeTerm<State>(t) {}
    FwdRTimeTerm(FwdRTimeTerm<State> const & t)
      : FwdTimeTerm<State>(t) {}
    ~FwdRTimeTerm() {}

    void setTargetTime(Time const & _t0) { 
      try {
	RealTime const & dt0=dynamic_cast<RealTime const &>(_t0);
	t0=dt0;
      }
      catch (bad_cast) {
	RVLException e;
	e<<"Error: FwdRTimeTerm::setTargetTime\n";
	throw e;
      }
    }
    void setTargetTime(double it) { t0=it; t0.write(cout); }
    Time const & getTargetTime() const { return t0; }

  };

  /** BwdTimeTerm fully implemented with StdDiscreteTime data member */
  template<typename State>
  class BwdRTimeTerm: public BwdTimeTerm<State> {
  private:
    RealTime t0;
    BwdRTimeTerm();
  public:
    BwdRTimeTerm(State const & t)
      : BwdTimeTerm<State>(t) {}
    BwdRTimeTerm(BwdRTimeTerm<State> const & t)
      : BwdTimeTerm<State>(t) {}
    ~BwdRTimeTerm() {}

    void setTargetTime(Time const & _t0) { 
      try {
	RealTime const & dt0=dynamic_cast<RealTime const &>(_t0);
	t0=dt0;
      }
      catch (bad_cast) {
	RVLException e;
	e<<"Error: BwdRTimeTerm::setTargetTime\n";
	throw e;
      }
    }
    void setTargetTime(double it) { t0=it; }
    Time const & getTargetTime() const { return t0; }

  };

}
#endif
