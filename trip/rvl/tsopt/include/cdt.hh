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
#ifndef __RVL_CDT
#define __RVL_CDT

#include "t.hh"

namespace TSOpt {

  using RVL::Writeable;
  using RVL::RVLException;

  /** Base time class for all fixed-step simulations */
  class DiscreteTime: public Time {
  public:
    DiscreteTime() {}
		  
    DiscreteTime(DiscreteTime const &) {}
    ~DiscreteTime() {  }

	/** Assignment - from int */
    using TSOpt::Time::operator=;
    virtual DiscreteTime & operator=(int it) = 0;

    /** increment operator - increments underlying int */
    virtual DiscreteTime & operator++() = 0;
    /** decrement operator - increments underlying int */
    virtual DiscreteTime & operator--() = 0;
   
    /** Assignment - inherited */
   // virtual DiscreteTime & operator=(DiscreteTime const & t) = 0;

              

  };


  /** StdDiscreteTime object - owns an int. Must be initialized by
      assignment. An overload is provided permitting direct assignment
      of the int data member. */

  class StdDiscreteTime: public DiscreteTime {
  private: 
    int * _it;
  public:
    StdDiscreteTime()
      : _it(NULL) {}
    StdDiscreteTime(StdDiscreteTime const & t) 
      : _it(NULL) { if (t._it) { _it=new int; *_it=*(t._it); }
    }
    ~StdDiscreteTime() { if (_it) { delete _it; _it=NULL; } }

    /** Assignment - inherited */
    Time & operator=(Time const & t);

    /** Assignment - avoid default */
    StdDiscreteTime & operator=(StdDiscreteTime const & t);

    /** Assignment - from int */
    StdDiscreteTime & operator=(int it);

    /** increment operator - increments underlying int */
    StdDiscreteTime & operator++();
 
	/** deccrement operator - increments underlying int */
    StdDiscreteTime & operator--();

    /** Comparison of int data */
    bool operator<(Time const &) const;
    /** Comparison of int data */
    bool operator>(Time const &) const;
    ostream & write(ostream & str) const;

    /** return integer time - by value, hence read-only */
    int getint() const { return *_it; }

	StdDiscreteTime * clone() const {
		return new StdDiscreteTime(*this);
	}
		  
    
  };

}

#endif
