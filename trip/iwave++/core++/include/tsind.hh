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
#ifndef __RVL_TSIND
#define __RVL_TSIND

#include "seamx_headers.hh"
#include "cdt.hh"

namespace TSOpt {

  using RVL::Writeable;
  using RVL::RVLException;

  class TSIndex: public DiscreteTime{
  private:
    // TSIndex();
  protected:
    mutable TIMESTEPINDEX lts;  /* local copy of time enables to make TSIndex array */
    TIMESTEPINDEX & ts;  /* need to associate with the time in simulation */
   /* ts and lts may hold different value, but all comparisons and operations
      with TSIndex based on ts */
                         
  public:
    TSIndex(): ts(lts){lts.it = 0; lts.iv = -1; lts.dt = 0;}
    TSIndex(TIMESTEPINDEX & _ts): ts(_ts) {
      lts.it = _ts.it; 
      lts.iv = _ts.iv;
      lts.dt = _ts.dt;
    }
    TSIndex(TSIndex const & t): ts(lts) {
      lts.it = t.ts.it;
      lts.iv = t.ts.iv; 
      lts.dt = t.ts.dt;
    }
    TSIndex * clone() const { return new TSIndex(*this); }
    ~TSIndex() {}
    /** assignment - post-construction initialization */
    Time & operator=(Time const &t);
    TSIndex & operator=(TSIndex const & t);
    TSIndex & operator=(int it);
    bool operator<(Time const & t) const;
    bool operator>(Time const & t) const;
    TSIndex & operator++();
    TSIndex & operator--();
    void synclts() {
      lts.it = ts.it; 
      lts.iv = ts.iv; 
      lts.dt = ts.dt;
    }
    // int getint() const { return getCstruct().it;} /* to make CPsim compiling*/
    ireal getTime() { return get_time(this->ts); }

    TIMESTEPINDEX & getCstruct() { return ts; }
    TIMESTEPINDEX const & getCstruct() const { return ts; }

    ostream & write(ostream & str) const;
  };

}

#endif
