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

#ifndef __RVL_LDCEVAL
#define __RVL_LDCEVAL

#include "localdata.hh"

#if __cplusplus >= 201103L   // We are using C++11 or a later version
#define NOEXCEPT_FALSE noexcept(false)
#else
#define NOEXCEPT_FALSE
#endif

namespace RVL {

  /** Evaluation type, mixin for definition of LFOs */
  template<class DataType>
  class LocalEvaluation {
  public:
    LocalEvaluation() {}
    LocalEvaluation(const LocalEvaluation<DataType> &) {}
    virtual ~LocalEvaluation() NOEXCEPT_FALSE {}
    /** Eval method for LDCs */    
    virtual void operator()(LocalDataContainer<DataType> & target,
			    vector<LocalDataContainer<DataType> const *> & sources) = 0;
  };

  /** Generic local function object, i.e. one that requires local
      data. One of many many possible implementations.
  */
  
  template<class DataType> 
  class LocalFunctionObject
    : public FunctionObject, public LocalEvaluation<DataType> {
  public:
    LocalFunctionObject() {}
    LocalFunctionObject(const LocalFunctionObject<DataType> &) {}
    virtual ~LocalFunctionObject() {}
  };
    
  /** Unary local evaluation mixin. */

  template<class DataType>
  class UnaryLocalEvaluation: public LocalEvaluation<DataType> {
  public:
    UnaryLocalEvaluation() {}
    UnaryLocalEvaluation(const UnaryLocalEvaluation<DataType> &) {}
    virtual ~UnaryLocalEvaluation() NOEXCEPT_FALSE {}

    /** Evaluation method for LDCs */
    virtual void operator () (LocalDataContainer<DataType> & target) = 0;

    /** Generic evaluation method */
    using LocalEvaluation<DataType>::operator();
    virtual void operator()(LocalDataContainer<DataType> & target,
			    vector<LocalDataContainer<DataType> const *> & sources) {
      try {
	if (sources.size() != 0) {
	  RVLException e;
	  e<<"Error: UnaryLocalEvaluation::operator() (generic)\n";
	  e<<"vector of sources not of length zero\n";
	  throw e;
	}
	(*this)(target);
      }
      catch (RVLException & e) {
	e<<"\ncalled from UnaryLocalEvaluation::operator() (generic)\n";
	throw e;
      }
    }
  };

  /** Unary local function object: takes one local data container
      argument, no return value. */

  template<class DataType>
  class UnaryLocalFunctionObject
    : public FunctionObject, public UnaryLocalEvaluation<DataType> {
  public:
    UnaryLocalFunctionObject() {}
    UnaryLocalFunctionObject(const UnaryLocalFunctionObject<DataType> &) {}
    virtual ~UnaryLocalFunctionObject() NOEXCEPT_FALSE {}
  };

  /** Binary local evaluation mixin. */

  template<class DataType>
  class BinaryLocalEvaluation: public LocalEvaluation<DataType> {
  public:
    BinaryLocalEvaluation() {}
    BinaryLocalEvaluation(const BinaryLocalEvaluation<DataType> &) {}
    virtual ~BinaryLocalEvaluation() NOEXCEPT_FALSE {}
    /** Evaluation method for LDCs */
    using LocalEvaluation<DataType>::operator();
    virtual void operator () (LocalDataContainer<DataType> & target,
			      LocalDataContainer<DataType> const & source) = 0;

    /** Generic evaluation method */
    void operator()(LocalDataContainer<DataType> & target,
		    vector<LocalDataContainer<DataType> const *> & sources) {
      try {
	if (sources.size() != 1) {
	  RVLException e;
	  e<<"Error: BinaryLocalFunctionObject::operator() (generic)\n";
	  e<<"vector of sources not of length zero\n";
	  throw e;
	}
	(*this)(target,*(sources[0]));
      }
      catch (RVLException & e) {
	e<<"\ncalled from BinaryLocalFunctionObject::operator() (generic)\n";
	throw e;
      }
    }
  };

  /** Binary local function object: takes two local data container
      arguments, no return value. */

  template<class DataType>
  class BinaryLocalFunctionObject
    : public FunctionObject, public BinaryLocalEvaluation<DataType> {
  public:
    BinaryLocalFunctionObject() {}
    BinaryLocalFunctionObject(const BinaryLocalFunctionObject<DataType> &) {}
    virtual ~BinaryLocalFunctionObject() NOEXCEPT_FALSE {}
  };

  /** Ternary local evaluation mixin. */

  template<class DataType>
  class TernaryLocalEvaluation: public LocalEvaluation<DataType> {
  public:
    TernaryLocalEvaluation() {}
    TernaryLocalEvaluation(const TernaryLocalEvaluation<DataType> &) {}
    virtual ~TernaryLocalEvaluation() {}

    /** Evaluation method for LDCs */
    using LocalEvaluation<DataType>::operator();
    virtual void operator () (LocalDataContainer<DataType> & target,
			      LocalDataContainer<DataType> const & source1,
			      LocalDataContainer<DataType> const & source2) = 0;

    /** Generic evaluation method */
    virtual void operator()(LocalDataContainer<DataType> & target,
			    vector<LocalDataContainer<DataType> const *> & sources) {
      try {
	if (sources.size() != 2) {
	  RVLException e;
	  e<<"Error: TernaryLocalFunctionObject::operator() (generic)\n";
	  e<<"vector of sources not of length 2\n";
	  throw e;
	}
	(*this)(target,*(sources[0]),*(sources[1]));
      }
      catch (RVLException & e) {
	e<<"\ncalled from TernaryLocalFunctionObject::operator() (generic)\n";
	throw e;
      }
    }
  };

  /** Ternary local function object: takes three local data container
      arguments, no return value. */

  template<class DataType>
  class TernaryLocalFunctionObject: 
    public FunctionObject, public TernaryLocalEvaluation<DataType> {
  public:
    TernaryLocalFunctionObject() {}
    TernaryLocalFunctionObject(const TernaryLocalFunctionObject<DataType> &) {}
    virtual ~TernaryLocalFunctionObject() {}
  };

  /** Quaternary local evaluation mixin. */

  template<class DataType>
  class QuaternaryLocalEvaluation: public LocalEvaluation<DataType> {
  public:
    QuaternaryLocalEvaluation() {}
    QuaternaryLocalEvaluation(const QuaternaryLocalEvaluation<DataType> &) {}
    virtual ~QuaternaryLocalEvaluation() {}

    /** Evaluation method for LDCs */
    using LocalEvaluation<DataType>::operator();
    virtual void operator () (LocalDataContainer<DataType> & target,
			      LocalDataContainer<DataType> const & source1,
			      LocalDataContainer<DataType> const & source2,
			      LocalDataContainer<DataType> const & source3) = 0;

    /** Generic evaluation method */
    virtual void operator()(LocalDataContainer<DataType> & target,
			    vector<LocalDataContainer<DataType> const *> & sources) {
      try {
	if (sources.size() != 3) {
	  RVLException e;
	  e<<"Error: QuaternaryLocalFunctionObject::operator() (generic)\n";
	  e<<"vector of sources not of length 3\n";
	  throw e;
	}
	(*this)(target,*(sources[0]),*(sources[1]),*(sources[2]));
      }
      catch (RVLException & e) {
	e<<"\ncalled from QuaternaryLocalFunctionObject::operator() (generic)\n";
	throw e;
      }
    }

  };

  /** Quaternary local function object: takes four local data container
      arguments, no return value. */

  template<class DataType>
  class QuaternaryLocalFunctionObject: 
    public FunctionObject, public QuaternaryLocalEvaluation<DataType> {
  public:
    QuaternaryLocalFunctionObject() {}
    QuaternaryLocalFunctionObject(const QuaternaryLocalFunctionObject<DataType> &) {}
    virtual ~QuaternaryLocalFunctionObject() {}
  };

}

#endif
