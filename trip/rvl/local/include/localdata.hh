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

#ifndef __RVL_LDC
#define __RVL_LDC

#include "data.hh"

namespace RVL {

/** LocalDataContainers provide direct access to data samples.
    Principal methods return the number of data samples (
    getSize()) and a pointer to Scalar (getData()). These and the
    reporting methods (two versions of write(...)) are pure virtual.

    LocalDataContainer abstracts the function of the STL valarray
    class, and is semantically identical (though of course
    syntactically incompatible) with many "vector" or "array"
    classes found in other OO numerics libraries.
    
    Since data access is uniform througout all LocalDataContainer
    child classes, it becomes possible to provide standard
    implementations of function objects performing common tasks, such
    as the basic linear algebra subroutines. A collection of these FOs
    form the RVLTools FO library (see functions.H) [this should be a
    link].

    The DataContainer eval(...) method is implemented in the base
    class DataContainer, via delegation to the operator() method of
    the various FunctionObject base class. For LocalFunctionObjects
    (defined below), this will involve downcasting all of the
    DataContainers which appear in the argument list to
    LocalDataContainers, as well as checking any restrictions on
    numbers of arguments. Thus LocalFunctionObjects can only be
    evaluated by LocalDataContainers. 

    It is also useful to provide a way to view a contiguous subarray
    of the data of a LocalDataContainer as a LocalDataContainer. This
    is the role of the LocalDataContainerSection class.

*/

  template<class DataType>
  class LocalEvaluation;

  template<class DataType>
  class LocalConstEval;

  template<class DataType>
  class LocalDataContainer: public DataContainer {

  public:

    LocalDataContainer() {}
    LocalDataContainer(const LocalDataContainer<DataType> & D) {}
    virtual ~LocalDataContainer() {}
  
    /** return size of local data container */
    virtual size_t getSize() const = 0;

    /** return address of writable data array */
    virtual DataType * getData() = 0;

    /** return address of read-only data array */
    virtual DataType const * getData() const = 0;

    /** local evaluation: defined at this level so that subtypes do not
	need to re-implement. The natural "stupid" implementation is
	the right one. 
    */
    void eval(FunctionObject & f,
	      vector<DataContainer const *> & x) {
      try {
	vector<LocalDataContainer<DataType> const *> lx(x.size());
	for (size_t i=0;i<x.size();i++) {
	  if (!(lx[i] = dynamic_cast<LocalDataContainer<DataType> const *>(x[i]))) {
	    RVLException e;
	    e<<"Error: LocalDataContainer::eval(FO,...)\n";
	    e<<"at least one of the input DCs is not an LDC\n";
	    throw e;
	  }
	}
	LocalEvaluation<DataType> & lf 
	  = dynamic_cast<LocalEvaluation<DataType> &>(f);
	lf(*this,lx);
      }
      catch (bad_cast) {
	RVLException e;
	e<<"Error: LocalDataContainer::eval(FO,...)\n";
	e<<"FO is not an LFO\n";
	throw e;
      }	
      catch (RVLException & e) {
	e<<"\ncalled from LocalDataContainer::eval\n";
	throw e;
      }
    }

    /** Similar evaluation method for FOCEs. */
    void eval(FunctionObjectConstEval & f,
	      vector<DataContainer const *> & x) const {
      try {
	vector<LocalDataContainer<DataType> const *> ex(x.size()+1);
	ex[0]=this;
	for (size_t i=0;i<x.size();i++) {
	  if (!(ex[i+1] = dynamic_cast<LocalDataContainer<DataType> const *>(x[i]))) {
	    RVLException e;
	    e<<"Error: LocalDataContainer::eval(FO,...)\n";
	    e<<"at least one of the input DCs is not an LDC\n";
	    throw e;
	  }
	}
	LocalConstEval<DataType> & lf 
	  = dynamic_cast<LocalConstEval<DataType> &>(f);
	lf(ex);
      }
      catch (RVLException & e) {
	e<<"\ncalled from LocalDataContainer::eval(...) const\n";
	throw e;
      }
    }
  };

  /** A factory class for LocalDataContainers. Useful in defining
      the partitioned data container class, also LocalSpaces. */

  template<class DataType>
  class LocalDataContainerFactory: public DataContainerFactory {

  public:

    LocalDataContainerFactory() {}
    LocalDataContainerFactory( LocalDataContainerFactory<DataType> &) {}
    virtual ~LocalDataContainerFactory() {}

    /** returns dynamically allocated LocalDataContainer */
    virtual LocalDataContainer<DataType> * buildLocal() const = 0;

    DataContainer * build() const { return buildLocal(); }

    /** returns size of LocalDataContainer product */
    virtual size_t getSize() const = 0;

  };


  /** A view of a contiguous portion of another LDC. */
  template<class DataType>
  class LocalDataContainerSection: public LocalDataContainer<DataType> {

  private:
    LocalDataContainer<DataType> & src;
    size_t begin;
    size_t length;

    LocalDataContainerSection() {}

  public:
    /** Copy constructor copies reference. Refers to same source LDC. */
    LocalDataContainerSection(const LocalDataContainerSection<DataType> & D)
      :src(D.src), begin(D.begin), length(D.length) {}
    /** Main constructor: initializes source LDC reference, begin index,
	length. Checks that section is legal. */
    LocalDataContainerSection(LocalDataContainer<DataType> & _src,
			      size_t _begin,
			      size_t _length)
      : src(_src), begin(_begin), length(_length) {
      if (begin+length>src.getSize()) {
	RVLException e;
	e<<"Error: LocalDataContainerSection constructor\n";
	e<<"section overruns end of source LDC\n";
	e<<"length of source LDC = "<<src.getSize()<<"\n";
	e<<"begin index = "<<begin<<" length of section = "<<length<<"\n";
	throw e;
      }	
    }
    virtual ~LocalDataContainerSection() {}
  
    size_t getSize() const { return length; }
    DataType * getData() { return &((src.getData())[begin]); }
    DataType const * getData() const { return &((src.getData())[begin]); }

    ostream & write(ostream & str) const {
      str<<"LocalDataContainerSection object\n";
      str<<"LDC view of section beginning at "<<begin<<" of length "<<length<<" of\n";
      str<<"elements:\n";
      for (size_t i=0;i<length;i++) {
	str<<"index = "<<i<<" value = "<<getData()[i]<<"\n";
      }
      return str;
    }
  };
  
}

#endif
