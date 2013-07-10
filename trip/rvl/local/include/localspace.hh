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

#ifndef __RVL_LOCSP
#define __RVL_LOCSP

#include "localdata.hh"
#include "space.hh"

namespace RVL {

  /** Abstract Space class for local data structures, i.e. data
      structures which can sensibly return access to their data via
      pointer to Scalar. Convenient to implement as StdSpace subclass.  */

  template<class Scalar, class DataType = Scalar>
  class LocalSpace: public StdSpace<Scalar, DataType> {

  protected:

    virtual  LocalDataContainerFactory<DataType> & getLDCF() const = 0;
    DataContainerFactory & getDCF() const { return getLDCF(); }

  public:

    LocalSpace() {}
    LocalSpace(const LocalSpace<Scalar, DataType> & sp) {}
    virtual ~LocalSpace() {}

    /** returns dynamically allocated LocalDataContainer */
    LocalDataContainer<DataType> * buildLocalDataContainer() const {
      return getLDCF().buildLocal();
    }

    DataContainer * buildDataContainer() const { 
      return buildLocalDataContainer(); 
    }

    bool isCompatible(DataContainer const & dc) const {
      try {
	return getLDCF().isCompatible(dc);
      }
      catch (RVLException & e) {
	e<<"\ncalled from RnSpace::isCompatible\n";
	throw e;
      }
    }
  };

  /** Vector in a LocalSpace, which can sensibly return a pointer
      to its data array - and does so (getData), also has access 
      to its with its dimension (getSize). 
  */

  template<class Scalar, class DataType = Scalar>
  class LocalVector: public Vector<Scalar> {

  private:

    LocalDataContainer<DataType> * ldc;

    LocalVector() {}

  public:

    /** Copy constructor. Uses Vector copy constructor, provides an
	internal DC pointer which is explicitly local. */
    LocalVector( const Vector<Scalar> & v )
      : Vector<Scalar>(v) {
      if (!(ldc = 
	    dynamic_cast<LocalDataContainer<DataType> *>(Vector<Scalar>::getDataContainer()))) {
	RVLException e;
	e<<"Error: LocalVector copy constructor\n";
	e<<"input not a local vector somehow\n";
	throw e;
      }
    }

    /** Constructor from Space. Provides an
	internal DC pointer which is explicitly local. */
    LocalVector(const Space<Scalar> & sp )
      : Vector<Scalar>(sp) {
      if (!(ldc = 
	    dynamic_cast<LocalDataContainer<DataType> *>(Vector<Scalar>::getDataContainer()))) {
	RVLException e;
	e<<"Error: LocalVector constructor\n";
	e<<"input space does not generate LocalDataContainers\n";
	throw e;
      }
    }

    ~LocalVector() {}

    /** returns dimension */
    size_t getSize() const { return ldc->getSize(); }
    /** returns data pointer (array) WWS 26.10.11: nonconst version
	increments version ref -- as it should, since implicitly
	invoking this method opens the internal state of the vector to
	change, as the nonconst return is an lvalue. IMPORTANT NOTE:
	C++ specifies that the least restrictive version, in this case
	non-const, of each method will be called, consistent with
	other conditions. Therefore almost any access to data will
	cause version increment, unless called on a const reference to
	the vector or within a const function which implicitly makes
	its reference to this type const.
     */
    DataType * getData() { this->Vector<Scalar>::incrementVersion(); return ldc->getData(); }
    DataType const * getData() const { return ldc->getData(); }

  };
}
#endif
