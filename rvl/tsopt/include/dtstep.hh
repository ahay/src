/*************************************************************************

Copyright Rice University, 2004, 2005, 2006, 2007
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
#ifndef __RVL_DTSIM
#define __RVL_DTSIM

#include "cdt.hh"
#include "sim.hh"

namespace TSOpt {

  /** Standard behaviour for next discrete time, forward. Implemented
      via StdDiscreteTime data member and increment operator.

      run method still abstract.
  */
  template<typename State>
  class StdDFTStep: public StdTimeStep<State> {

  private:
    
    mutable StdDiscreteTime tnext;

  public:

    StdDFTStep(): StdTimeStep<State>() {}
    StdDFTStep(StdDFTStep<State> const & t)
      : StdTimeStep<State>(t), tnext(t.tnext) {}
    virtual ~StdDFTStep() {}

    Time const & getNextTime() const {
      tnext=this->getTime();
      ++tnext;
      return tnext;
    }
    
  };

  /** Standard behaviour for next discrete time, backward. Implemented
      via StdDiscreteTime data member and deccrement operator.

      run method still abstract.
  */
  template<typename State>
  class StdDBTStep: public StdTimeStep<State> {

    mutable StdDiscreteTime tnext;

  public:

    StdDBTStep(): StdTimeStep<State>() {}
    StdDBTStep(StdDBTStep<State> const & t)
      : StdTimeStep<State>(t), tnext(t.tnext) {}
    virtual ~StdDBTStep() {}

    Time const & getNextTime() const {
      tnext=this->getTime();
      --tnext;
      return tnext;
    }
  };
}

#endif
