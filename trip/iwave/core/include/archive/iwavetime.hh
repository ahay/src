/*************************************************************************

Copyright Rice University, 2004-2014
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
#ifndef __IWAVETIME
#define __IWAVETIME

#include "seamx_headers.hh"
#include "cdt.hh"
#include "tterm.hh"

namespace TSOpt {

  using RVL::Writeable;
  using RVL::RVLException;

  class IWaveTime: public DiscreteTime{
  private:
    IWaveTime();
  protected:

    TIMESTEPINDEX & ts;  /* need to associate with the time in simulation */
   /* ts and lts may hold different value, but all comparisons and operations
      with IWaveTime based on ts */
                         
  public:

    IWaveTime(TIMESTEPINDEX & _ts): ts(_ts) {}
    IWaveTime(IWaveTime const & t): ts(t.ts) {}
    IWaveTime * clone() const { return new IWaveTime(*this); }
    ~IWaveTime() {}
    /** assignment - post-construction initialization */
    Time & operator=(Time const &t);
    IWaveTime & operator=(IWaveTime const & t);
    IWaveTime & operator=(int it);
    bool operator<(Time const & t) const;
    bool operator>(Time const & t) const;
    IWaveTime & operator++();
    IWaveTime & operator--();

    // int getint() const { return getCstruct().it;} /* to make CPsim compiling*/
    ireal getTime() { return get_time(this->ts); }

    TIMESTEPINDEX & getCstruct() { return ts; }
    TIMESTEPINDEX const & getCstruct() const { return ts; }

    ostream & write(ostream & str) const;
  };
 
  class IWaveTimeTerm: public TimeTerm {
  private:

    bool fwd;
    mutable IWaveTime t0;
    IWaveTime const & t;
    IWaveTimeTerm();
    IWaveTimeTerm(IWaveTimeTerm const &);

  public:

    IWaveTimeTerm(IWaveTime const & _t,
		  IWaveTime _t0,
		  bool _fwd = true) 
      : t(_t), t0(_t0), fwd(_fwd) {}
    void setTargetTime(Time const & _t0) { t0 = _t0; }
    Time const & getTargetTime() const { return t0; }
    bool query() {
      if (fwd) return ((t<t0) || (t==t0));
      else return ((t>t0) || (t==t0));
    }
    ostream & write(ostream & str) const {
      str<<"IWaveTimeTerm\n";
      str<<"current time = \n";
      t.write(str);
      str<<"target time  = \n";
      t0.write(str);
      return str;
    }
  };
}

#endif
