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
#ifndef __RVL_POLICY
#define __RVL_POLICY

#include "except.hh"

namespace TSOpt {

  using RVL::RVLException;

  /** generic creation policy using new, after Alexandrescu.
      Only nontrivial method is create. Types:
      <ul>
      <li> T: type to be created; </li>
      <li> C: data type for constructor </li>
      </ul>
      Of course for this to work there must be a constructor for
      type T that takes a C; this may require adding a struct and 
      a new constructor to an existing class.
  */

  template<typename C, typename T>
  class OpNewCreatePolicy {
  protected:
    T * create(C & c) const {
      try {
	return new T(c);
      }
      catch (RVLException & e) {
	e<<"\ncalled from opNewCreatePolicy::create\n";
	throw e;
      }
    }
  };

  template<typename C, typename T>
  class OpNewConstCreatePolicy {
  protected:
    T * create(C const & c) const {
      try {
	return new T(c);
      }
      catch (RVLException & e) {
	e<<"\ncalled from OpNewConstCreatePolicy::create\n";
	throw e;
      }
    }
  };

  /* for test purposes - nada */  
  template<typename C, typename T>
  class PolicyBase {
  public:
      virtual ~PolicyBase() {}
  protected:
    virtual T * create(C & c) const {
      RVLException e;
      e<<"Error: PolicyBase::create\n";
      e<<"it's a really bad idea to actually call this method\n";
      e<<"you need to override it in a subclass\n";
      throw e;
    }
  };

  template<typename C, typename T>
  class ConstPolicyBase {
  protected:
    virtual T * create(C const & c) const {
      RVLException e;
		e<<"Error: ConstPolicyBase::create\n";
      e<<"it's a really bad idea to actually call this method\n";
      e<<"you need to override it in a subclass\n";
      throw e;
    }
  };

  /** also for test purposes - allows public access to policy
      create methods */

  template<typename P, typename C, typename T>
  class Creator: public P {
  public:
    T * create(C & c) const {
      try {
	return P::create(c); 
      }
      catch (RVLException & e) {
	e<<"called from Creator::creat\n";
	throw e;
      }
    }
  };

  template<typename P, typename C, typename T>
  class ConstCreator: public P {
  public:
    T * create(C const & c) const {
      try {
	return P::create(c); 
      }
      catch (RVLException & e) {
	e<<"called from ConstCreator::creat\n";
	throw e;
      }
    }
  };
}

#endif
