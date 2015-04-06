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

#ifndef __RVLALG_UMIN_LBFGSDS_H
#define __RVLALG_UMIN_LBFGSDS_H

#include "umintable.hh"
#include "lbfgsalg.hh"
#include "lnsrchDS.hh"

namespace RVLUmin{
  using namespace RVLAlg;
  using RVL::Vector;
  using RVL::Functional;
  using RVL::FunctionalEvaluation;

  template<class Scalar>
  class UMinLBFGS_DS: public StateAlg< Vector<Scalar> > {

  private:

    LineSearchDS<Scalar> ds;
    UMinLBFGS<Scalar> umin;
    bool display;
    
    ostream & str;

    UMinLBFGS_DS();

  public:

    UMinLBFGS_DS(Functional<Scalar> & f,
		 Vector<Scalar> & x,
		 UMinTable<Scalar> & holder,
		 ostream & _str=cout)
      : ds(x.getSpace(),holder,_str), 
	umin(f,x,ds,holder,_str),
	str(_str) {
      int val = 0;
      if (holder.getValue("DispFlag",val)) display=false;
      if (val) display=true;
    }

    UMinLBFGS_DS(const UMinLBFGS_DS<Scalar> & bu)
      : ds(bu.ds), umin(bu.umin), str(bu.str), display(bu.display) {}

    virtual ~UMinLBFGS_DS() {}

    void setState(const Vector<Scalar> & s) { umin.setState(s); }

    Vector<Scalar> & getState() { return umin.getState(); }

    bool run() {
      try {
	if (display) {
	  str<<"/***************************************************\n";
	  str<<" *   Unconstrained Minimisation by Limited Memory  *\n";
	  str<<" * BFGS Iteration with Dennis-Schnabel Line Search *\n";
	  str<<" ***************************************************/\n";
	}
	return umin.run();
      }
      catch (RVLException & e) {
	e<<"\ncalled from UMinLBFGS_BT\n";
	throw e;
      }
    }
  };
}

#endif
