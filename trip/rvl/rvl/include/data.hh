/*************************************************************************

Copyright Rice University, 2004, 2005, 2006.
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

#ifndef __RVL_DC
#define __RVL_DC

#include "utility.hh"

namespace RVL {

  /* Main interfaces */

  /** The FunctionObject / DataContainer hiearchy follows the Acyclic
      Visitor design pattern (see Martin 2002 for more on this). The
      base visitor class FunctionObject has no (nontrivial) methods,
      which breaks the dependency between visitors and visiteds
      inherent in the GoF Visitor pattern. The meaning of such an
      object is entirely defined by its use, for which see
      DataContainer, below, and Vector in Space.H.
      
      Concrete subtypes of FunctionObject will typically apply runtime
      tests to constrain the types of data container source and target
      arguments, as well as the number of arguments.

  */

  class FunctionObject: public Writeable {

  public:
    
    FunctionObject() {}
    FunctionObject(const FunctionObject &) {}
    virtual ~FunctionObject() {}

    /** Name method - returns string. Used in standard implementations
	of reporting methods, which simply report the name of the
	function object. */
    virtual string getName() const = 0;
    
    /** report to ostream - can be overridden in subclasses. */
    virtual ostream & write(ostream & str) const {
      str<<"Function Object"<<"\n";
      str<<"  name = "<<getName()<<"\n";
      return str;
    }

  };

  /** Membership
  */

  class FunctionObjectConstEval: public Writeable {
  public: 
    virtual ~FunctionObjectConstEval() {}
    virtual string getName() const = 0;
    virtual ostream & write(ostream & str) const {
      str<<"Function Object for Constant Evaluation:"<<"\n";
      str<<"  name = "<<getName()<<"\n";
      return str;
    }
  };

  /** Mixin class for types that store, update, and return "small"
      objects, meaning effectively those with usable copy
      semantics. Appends a simple implemented container attribute,
      which has the added effect of testing for the existence of
      assignment at compile time. Intended for use with FOCE to create
      reduction function objects.
  */

  template<typename Scalar>
  class ScalarRedn {
  private:
    mutable Scalar val;
    /** since there is no universal default value of type Scalar, there is
	no natural definition of a default constructor for this class */
    ScalarRedn();
  public:
    /** main constructor. Note that bytewise copy construction is exactly 
	what is intended, so let the compiler do it. */
    ScalarRedn(Scalar _val) 
      : val(_val) {}
#if __cplusplus >= 201103L                // We are using C++11 or a later version      
    virtual ~ScalarRedn() noexcept(false) {}
#else
    virtual ~ScalarRedn() {}
#endif

    /** post-construction (re)initialization - default is undefined,
	so pure virtual */
    virtual void setValue() = 0;
    /** post-construction (re)initialization */
    virtual void setValue(Scalar _val) { val=_val; }
    /** access - virtual so that additional behaviour may be added 
	in child class overrides, for expl cross-process reduction */
    virtual Scalar getValue() const { return val; }
  };
  
  /** Function object with const eval and scalar reduction attributes */
  template<typename ValType> 
  class FunctionObjectScalarRedn:
    public FunctionObjectConstEval, public ScalarRedn<ValType> {
  private:
    FunctionObjectScalarRedn();
  public:
    FunctionObjectScalarRedn(ValType val): ScalarRedn<ValType>(val) {}
    FunctionObjectScalarRedn(FunctionObjectScalarRedn<ValType> & f)
      : ScalarRedn<ValType>(f) {}
#if __cplusplus >= 201103L  // We are using C++11 or a later version
    virtual ~FunctionObjectScalarRedn() noexcept(false) {}
#else
    virtual ~FunctionObjectScalarRedn() {}
#endif
  };


  /** DataContainer is the principal RVL abstraction for types
      encapsulating data structures. This abstract type is intended to
      encompass an arbitrary number of levels of indirection (cartesian
      product structures, for example, or trees, or...) hence provides
      no actual access to data. All interaction with data takes place
      through evaluation of FunctionObjects.

      In GoF terms, the FO and DC hierarchies together exemplify the
      Visitor pattern, in its Acyclic variant defined by Martin. The
      DCs are the Visited types (i.e. the Elements, in terms of the
      GoF chapter).

      We have chosen to locate control of evaluation in the
      DataContainer class rather than in the FunctionObject
      class. This decision is standard in the Visitor pattern,
      i.e. that the visited type controls access, and allows
      FunctionObjects to be reused across wide categories of
      DataContainers, without needing to be aware of data structure
      details. If we were to locate control in the FunctionObject
      type, then each FO would need to be coded to deal correctly with
      every type of DC. Many FOs are simple scalar operations, like
      linear combination, which are evaluated in the same way on every
      DataContainer at the level of data, but which need not be aware
      of the overall structure if DataContainers control execution.

      Ultimately some DataContainer subtype must actually provide
      access to data, else nothing happens. An example of such a
      construction is the LocalDataContainer class template in
      LocalRVL.
  */

  class DataContainer: public Writeable {
  public:

    DataContainer() {}
    DataContainer(const DataContainer &) {}
    virtual ~DataContainer() {}

    /** Evaluate a function object.Concrete subtypes must supply this
	definition in order to properly exploit the data structures
	which they encapsulate. For example the ProductDataContainer
	structure defined in productdata.H implements evaluation via a
	loop over factors.
    */
    virtual void eval(FunctionObject & f, 
		      vector<DataContainer const *> & x) = 0;
	
    /** Const version of evaluation */
    virtual void eval(FunctionObjectConstEval & f, 
		      vector<DataContainer const *> & x) const = 0;
  };

  /** Factory class for DataContainers. Utmost simplicity. Typical
      use: implement Space::buildDataContainer in Space subclass. */

  class DataContainerFactory: 
    public Factory<DataContainer> {

  public:

    DataContainerFactory() {}
    DataContainerFactory(const DataContainerFactory &) {}
    virtual ~DataContainerFactory() {}
 
    /* returns dynamically allocated DataContainer */
    /* inherited from Factory - 25.04.10  
    virtual DataContainer * build() const = 0;
    */

    /** determines whether another DataContainerFactory builds the same
	kind of DataContainer. Usually this means effectively that the two
	are of the same type, determined by runtime type checking. Returns
	zero if different types, otherwise nonzero.
    */
    virtual bool compare( DataContainerFactory const & dcf) const = 0;

    /** determines whether a DataContainer is of the type built by this
	factory. Usually implemented through runtime type checking.
	Returns zero if not, else nonzero. 
    */
    virtual bool isCompatible(DataContainer const & dc) const = 0;

  };

}

#endif










