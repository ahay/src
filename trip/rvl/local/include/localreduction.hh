/*************************************************************************

Copyright Rice University, 2004-2011.
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

#ifndef __RVL_LDCREDN
#define __RVL_LDCREDN

#include "localdata.hh"

#if __cplusplus >= 201103L // We are using C++11 or a later version
#define NOEXCEPT_FALSE noexcept(false)
#else
#define NOEXCEPT_FALSE
#endif

namespace RVL {
  
  /** Local FOCE */
  template<typename DataType>
  class LocalConstEval {
  public:
    virtual ~LocalConstEval() NOEXCEPT_FALSE {}
    /** Eval method for LDCs */    
    virtual void operator()(vector<LocalDataContainer<DataType> const *> & sources) = 0;
  };

  /** Unary Local FOCE with templated return type interface.  */
  template<typename DataType>
  class UnaryLocalConstEval
    : public LocalConstEval<DataType> {
  public:
    virtual ~UnaryLocalConstEval() {}

    /** Evaluation method for LDCs */
    virtual void operator()(LocalDataContainer<DataType> const & source) = 0;

    /** Generic evaluation method */
    void operator()(vector<LocalDataContainer<DataType> const *> & sources) {
      try {
	if (sources.size() != 1) {
	  RVLException e;
	  e<<"Error: UnaryLocalConstEval::opeartor() (generic)\n";
	  e<<"input length != 1\n";
	  throw e;
	}
	return (*this)(*(sources[0]));
      }
      catch (RVLException & e) {
	e<<"\ncalled from UnaryLocalConstEval::operator() (generic)\n";
	throw e;
      }
    }

  };

  /** Binary Local FOCE with templated return type interface.  */
  template<typename DataType>
  class BinaryLocalConstEval
    : public LocalConstEval<DataType>  {
  public:
    virtual ~BinaryLocalConstEval() NOEXCEPT_FALSE {}

    /** Evaluation method for LDCs */
    virtual void operator () (LocalDataContainer<DataType> const & source1,
			      LocalDataContainer<DataType> const & source2) = 0;

    /** Generic evaluation method */
    void operator()(vector<LocalDataContainer<DataType> const *> & sources) {
      try {
	if (sources.size() != 2) {
	  RVLException e;
	  e<<"Error: BinaryLocalConstEval::opeartor() (generic)\n";
	  e<<"input length != 2\n";
	  throw e;
	}
	(*this)(*(sources[0]),*(sources[1]));
      }
      catch (RVLException & e) {
	e<<"\ncalled from BinaryLocalConstEval::operator() (generic)\n";
	throw e;
      }
    }
  };

  /** Ternary Local FOCE with templated return type interface.  */
  template<typename DataType, typename ValType = DataType>
  class TernaryLocalConstEval
    : public LocalConstEval<DataType> {
  public:
    virtual ~TernaryLocalConstEval() {}

    /** Evaluation method for LDCs */
    virtual void operator () (LocalDataContainer<DataType> const & source1,
			      LocalDataContainer<DataType> const & source2,
			      LocalDataContainer<DataType> const & source3) = 0;

    /** Generic evaluation method */
    void operator()(vector<LocalDataContainer<DataType> const *> & sources) {
      try {
	if (sources.size() != 3) {
	  RVLException e;
	  e<<"Error: TernaryLocalConstEval::operator() (generic)\n";
	  e<<"input length != 3\n";
	  throw e;
	}
	(*this)(*(sources[0]),*(sources[1]),*(sources[2]));
      }
      catch (RVLException & e) {
	e<<"\ncalled from TernaryLocalConstEval::operator() (generic)\n";
	throw e;
      }
    }

  };

  /** Quaternary Local FOCE with templated return type interface.  */
  template<typename DataType>
  class QuaternaryLocalConstEval
    : public LocalConstEval<DataType> {
  public:
    virtual ~QuaternaryLocalConstEval() {}

    /** Evaluation method for LDCs */
    virtual void operator () (LocalDataContainer<DataType> const & source1,
			      LocalDataContainer<DataType> const & source2,
			      LocalDataContainer<DataType> const & source3,
			      LocalDataContainer<DataType> const & source4) = 0;

    /** Generic evaluation method */
    void operator()(vector<LocalDataContainer<DataType> const *> & sources) {
      try {
	if (sources.size() != 4) {
	  RVLException e;
	  e<<"Error: QuaternaryLocalConstEval::opeartor() (generic)\n";
	  e<<"input length != 3\n";
	  throw e;
	}
	(*this)(*(sources[0]),*(sources[1]),*(sources[2]),*(sources[3]));
      }
      catch (RVLException & e) {
	e<<"\ncalled from QuaternaryLocalConstEval::operator() (generic)\n";
	throw e;
      }
    }
  };

  /** local const eval hierarchy. use for evaluations without
      templated reduction types - for example, when reduction target
      is external, accessed by reference. */

  template<typename DataType>
  class UnaryLocalFunctionObjectConstEval
    : public FunctionObjectConstEval, public UnaryLocalConstEval<DataType> {
  public:
    virtual ~UnaryLocalFunctionObjectConstEval() {}
  };

  template<typename DataType>
  class BinaryLocalFunctionObjectConstEval
    : public FunctionObjectConstEval, public BinaryLocalConstEval<DataType> {
  public:
    virtual ~BinaryLocalFunctionObjectConstEval() {}
  };

  template<typename DataType>
  class TernaryLocalFunctionObjectConstEval
    : public FunctionObjectConstEval, public TernaryLocalConstEval<DataType> {
  public:
    virtual ~TernaryLocalFunctionObjectConstEval() {}
  };

  template<typename DataType>
  class QuaternaryLocalFunctionObjectConstEval
    : public FunctionObjectConstEval, public QuaternaryLocalConstEval<DataType> {
  public:
    virtual ~QuaternaryLocalFunctionObjectConstEval() {}
  };

  /** local reduction function object hierarchy */

  template<typename DataType, typename ValType>
  class UnaryLocalFunctionObjectScalarRedn
    : public FunctionObjectScalarRedn<ValType>, public UnaryLocalConstEval<DataType> {
  private:
    UnaryLocalFunctionObjectScalarRedn();
  public:
    UnaryLocalFunctionObjectScalarRedn(ValType val)
      : FunctionObjectScalarRedn<ValType>(val) {}
    UnaryLocalFunctionObjectScalarRedn(UnaryLocalFunctionObjectScalarRedn<DataType,ValType> const & f) 
      : FunctionObjectScalarRedn<ValType>(f) {}
    virtual ~UnaryLocalFunctionObjectScalarRedn() {}
  };

  template<typename DataType, typename ValType>
  class BinaryLocalFunctionObjectScalarRedn
    : public FunctionObjectScalarRedn<ValType>, public BinaryLocalConstEval<DataType> {
  private:
    BinaryLocalFunctionObjectScalarRedn();
  public:
    BinaryLocalFunctionObjectScalarRedn(ValType val)
      : FunctionObjectScalarRedn<ValType>(val) {}
    BinaryLocalFunctionObjectScalarRedn(BinaryLocalFunctionObjectScalarRedn<DataType,ValType> const & f) 
      : FunctionObjectScalarRedn<ValType>(f) {}
    virtual ~BinaryLocalFunctionObjectScalarRedn() NOEXCEPT_FALSE {}
  };

  template<typename DataType, typename ValType>
  class TernaryLocalFunctionObjectScalarRedn
    : public FunctionObjectScalarRedn<ValType>, public TernaryLocalConstEval<DataType> {
  private:
    TernaryLocalFunctionObjectScalarRedn();
  public:
    TernaryLocalFunctionObjectScalarRedn(ValType val)
      : FunctionObjectScalarRedn<ValType>(val) {}
    TernaryLocalFunctionObjectScalarRedn(TernaryLocalFunctionObjectScalarRedn<DataType,ValType> const & f) 
      : FunctionObjectScalarRedn<ValType>(f) {}
    virtual ~TernaryLocalFunctionObjectScalarRedn() {}
  };

  template<typename DataType, typename ValType>
  class QuaternaryLocalFunctionObjectScalarRedn
    : public FunctionObjectScalarRedn<ValType>, public QuaternaryLocalConstEval<DataType> {
  private:
    QuaternaryLocalFunctionObjectScalarRedn();
  public:
    QuaternaryLocalFunctionObjectScalarRedn(ValType val)
      : FunctionObjectScalarRedn<ValType>(val) {}
    QuaternaryLocalFunctionObjectScalarRedn(QuaternaryLocalFunctionObjectScalarRedn<DataType,ValType> const & f) 
      : FunctionObjectScalarRedn<ValType>(f) {}
    virtual ~QuaternaryLocalFunctionObjectScalarRedn() {}
  };

}

#endif
