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
#ifndef __RVL_TTERM
#define __RVL_TTERM

#include "t.hh"
#include "terminator.hh"
#include "write.hh"

namespace TSOpt {

  using RVL::Writeable;
  using RVL::RVLException;
  using RVLAlg::Terminator;

  /** Time-aware terminator - facilities to set and get target time */
  class TimeTerm: public Terminator, public Writeable {
  public:
    TimeTerm() {}
    TimeTerm(TimeTerm const & t) {}
    virtual ~TimeTerm() {}

    /** the set method provides the only interface for
	alterning the Time object exposed by the get method */
    virtual void setTargetTime(Time const & _t0) = 0;
    /** get returns a read-only Time */
    virtual Time const & getTargetTime() const = 0;

  };

  /** monitors the Time encapsulated in a State, and returns
      true when the target Time managed by the TimeTerm is greater
      or equal. */
  template<typename State>
  class FwdTimeTerm: public TimeTerm {
  private:
    bool verbose;
    ostream & str;
  protected:
    State const & s;
  public:
    /** main constructor - setting the verbosity flag causes
	the query method to dump useful information to the 
	stream provided */
    FwdTimeTerm(State const & t,
		bool _verbose=false,
		ostream & _str=cout)
	: verbose(_verbose), str(_str), s(t) {}
    FwdTimeTerm(FwdTimeTerm<State> const & t)
	: verbose(t.verbose), str(t.str), s(t.s) {}
    virtual ~FwdTimeTerm() {}    
    virtual bool query() { 
      try {

	if (verbose) this->write(str);

	if (this->getTargetTime() <  s.State::getTime() ||
	    this->getTargetTime() == s.State::getTime()) {
	  return true; 
	}
	return false;
      }
      catch (RVLException & e) {
	e<<"\ncalled from FwdTimeTerm::query\n";
	throw e;
      }
    }  
    
    ostream & write(ostream & fp) const {
      fp<<"--------- beg FwdTimeTerm::write ---------------\n";
      fp<<"target time = ";
      this->getTargetTime().write(fp);
      fp<<"state time  = ";
      s.getTime().write(fp);
      fp<<"--------- end FwdTimeTerm::write ---------------\n";
      return fp;
    }

  };

  /** monitors the Time encapsulated in a State, and returns
      true when the target Time managed by the TimeTerm is lesser
      or equal. */
  template<typename State>
  class BwdTimeTerm: public TimeTerm {
  private:
    bool verbose;
    ostream & str;
  protected:
    State const & s;
  public:
    BwdTimeTerm(State const & t,
		bool _verbose=false,
		ostream & _str=cout)
	: verbose(_verbose), str(_str), s(t) {}
    BwdTimeTerm(BwdTimeTerm<State> const & t)
	: verbose(t.verbose), str(t.str), s(t.s) {}
    virtual ~BwdTimeTerm() {}    
    virtual bool query() { 
      try {
	if (verbose) this->write(str);
	if (this->getTargetTime() >  s.State::getTime() ||
	    this->getTargetTime() == s.State::getTime()) {
	  return true; 
	}
	return false;
      }
      catch (RVLException & e) {
	e<<"\ncalled from BwdTimeTerm::query\n";
	throw e;
      }
    }  
    
    ostream & write(ostream & fp) const {
      fp<<"--------- beg BwdTimeTerm::write ---------------\n";
      fp<<"target time = ";
      this->getTargetTime().write(fp);
      fp<<"state time  = ";
      s.getTime().write(fp);
      fp<<"--------- end BwdTimeTerm::write ---------------\n";
      return fp;
    }

  };

  /** Modeled after Alg::OrTerminator; gets TimeTerm services
      from first term */
  class OrTimeTerm: public TimeTerm {
    
  protected:
    TimeTerm & first;
    Terminator & second;
    OrTimeTerm();
  public:
    OrTimeTerm( TimeTerm &a, Terminator &b): first(a), second(b) {}
    OrTimeTerm( OrTimeTerm const & aa ): first(aa.first), second(aa.second) {}
    ~OrTimeTerm() {}
    virtual bool query() { return (first.query() || second.query()); }
    void setTargetTime(Time const & _t0) { first.setTargetTime(_t0); }
    Time const & getTargetTime() const { return first.getTargetTime(); }
    ostream & write(ostream & str) const; 
  };

  /** Modeled after Alg::AndTerminator; gets TimeTerm services
      from first term */
  class AndTimeTerm: public TimeTerm {
    
  protected:
    TimeTerm & first;
    Terminator & second;
    AndTimeTerm();
  public:
    AndTimeTerm( TimeTerm &a, Terminator &b): first(a), second(b) {}
    AndTimeTerm( AndTimeTerm const & aa ): first(aa.first), second(aa.second) {}
    ~AndTimeTerm() {}
    virtual bool query() { return (first.query() && second.query()); }
    void setTargetTime(Time const & _t0) { first.setTargetTime(_t0); }
    Time const & getTargetTime() const { return first.getTargetTime(); }
    ostream & write(ostream & str) const; 
  };

}
#endif
