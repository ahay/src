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

#ifndef __RVLALG_BOOL_TERMINATOR
#define __RVLALG_BOOL_TERMINATOR

#include "alg.hh"

/**  boolterm.H
     All Terminators which manipulate other terminators using boolean logic.

*/

namespace RVLAlg {

  /** A terminator which simply acts as a container for a boolean */
  class BoolTerminator: public Terminator {
  private:
    mutable bool ans;
  public:
    BoolTerminator(bool _ans=false): ans(_ans) {}
    BoolTerminator(BoolTerminator & bt): ans(bt.ans) {}
    ~BoolTerminator() {}
    void setValue() { ans=false; }
    void setValue(bool _ans) { ans=_ans; }
    bool query() { return ans; }
  };

  /** Build a new Terminator by combining the results of two others with 
      a logical AND.

      Currently obeys the short-circuit evaluation behavior in the C++ standard
      for logical operations.  Thus, if the first.query is false (continue looping)
      then the second query never occurs.

      Note that this object does NOT own the other two Terminators,
      and thus deleting one of them will cause nasty results.

      Ideally, it would be nice to have some sort of operator behavior
      so you could write

      term = (A AND B AND C) OR D

      where A,B,C,D were all Terminators, but I don't know how to do this yet.
      Of course, in most cases, one can get the desired results just by 
      combining the query calls with the proper logical operations, but I can
      forsee cases where we might need to pass a single Terminator into 
      an algorithm, in which case being able to build one out of combinations
      of others would be useful.
  */
  class AndTerminator: public Terminator {
  public:
    
    AndTerminator( Terminator &a, Terminator &b): first(a), second(b){}
    ~AndTerminator() {}
    
    virtual bool query() { 
      return (first.query() && second.query());
    }
    
  protected:
    Terminator & first;
    Terminator & second;
    
  };
  
  /** Build a new Terminator by combining the results of two others with 
      a logical OR.
      
      Currently obeys the short-circuit evaluation behavior in the C++ standard
      for logical operations. 
  */
  class OrTerminator: public Terminator {
  public:

    OrTerminator( Terminator &a, Terminator &b): first(a), second(b){}
    ~OrTerminator() {}

    virtual bool query() { 
      return (first.query() || second.query());
    }

  protected:
    Terminator & first;
    Terminator & second;
  };

  /** Build a new Terminator by inverting the result of another with 
      a logical NOT.
  */
  class NotTerminator: public Terminator {
  public:

    NotTerminator( Terminator &a): first(a){}
    ~NotTerminator() {}

    virtual bool query() { 
      return (!first.query()); 
    }

  protected:
    Terminator & first;
  };

  /** Build a new Terminator by combining the results of two others with 
      a logical XOR.

      Currently obeys the short-circuit evaluation behavior in the C++ standard
      for logical operations. 
  */
  class XorTerminator: public Terminator {
  public:

    XorTerminator( Terminator &a, Terminator &b): first(a), second(b){}
    ~XorTerminator() {}

    virtual bool query() {
      bool A = first.query();
      bool B = second.query();

      return ((A || B) && (!(A && B)));
    }

  protected:
    Terminator & first;
    Terminator & second;
  };

}

#endif
