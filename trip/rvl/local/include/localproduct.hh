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

#ifndef __RVL_PLDC
#define __RVL_PLDC


#include "localdata.hh"
#include "productdata.hh"

namespace RVL {

  /** Specialization of ProductDataContainer to case where all factors
      are local data containers. In designing such a class, it is imperative
      to decide what its primary identity shall be, in order to attach
      the class to the correct branch of the inheritance tree. Generally,
      the right decision is probably to attach it as "low down" as possible,
      i.e. with as many properties. I have followed that rule here: rather 
      than deriving this class from ProductDataContainer, which simply 
      implements virtual PDC properties, I have derived it from LocalDC,
      to which it adds properties. It's easy enough to write a wrapper 
      class to create a PDC out of one of these:
  */

  template<class DataType>
  class ProductLocalDataContainer: public LocalDataContainer<DataType> {

  public:

    ProductLocalDataContainer() {}
    ProductLocalDataContainer( ProductLocalDataContainer<DataType> & ) {}
    ~ProductLocalDataContainer() {}

    /** return number of components - note that this cannot be 
	getSize(), which means something else! */
    virtual int getNumberOfComponents() const = 0;
  
    /** return reference to ith component as LocalDataContainer */
    virtual LocalDataContainer<DataType> & operator[](int i) = 0;
    virtual LocalDataContainer<DataType> const & operator[](int i) const = 0;

  };

  /** wrapper class which turns a ProductLocalDataContainer
      into a ProductDataContainer */
  
  template<class DataType>
  class ProductDataContainerLDC: public ProductDataContainer {

  private:

    ProductLocalDataContainer<DataType> & pldc;

    ProductDataContainerLDC();

  public:

    ProductDataContainerLDC( ProductLocalDataContainer<DataType> & _pldc)
      : pldc(_pldc) {}
    ProductDataContainerLDC( const ProductDataContainerLDC<DataType> & pdldc )
      : pldc(pdldc.pldc) {}
    virtual ~ProductDataContainerLDC() {}

    virtual int getSize() const { return pldc.getNumberOfComponents(); }

    virtual DataContainer & operator[](int i) { 
      try {
	return pldc[i];
      }
      catch (RVLException & e) {
	e<<"\ncalled from ProductDataContainerLDC::operator[]\n";
	throw e;
      }
    }
    virtual DataContainer const & operator[](int i) const { 
      try {
	return pldc[i];
      }
      catch (RVLException & e) {
	e<<"\ncalled from ProductDataContainerLDC::operator[]\n";
	throw e;
      }
    }
  };

  /** ProductLocalDataContainer made by dicing up a LocalDataContainer.
      It is possible to treat the data as either a
      single array or as a block structured array. */

  template<class DataType>
  class PartitionedLocalDataContainer
    : public ProductLocalDataContainer<DataType> {

  private:

    LocalDataContainer<DataType> & ldc;
    vector<int> k;

    vector<LocalDataContainerSection<DataType> *> comps;

    PartitionedLocalDataContainer() {}

  public:

    /** Copy constructor. Makes new vector of pointers to the components. */
    PartitionedLocalDataContainer(const PartitionedLocalDataContainer<DataType> & p)
      : ldc(p.ldc), k(p.k), comps(k.size()) {
      try {
	int nblk=k.size();
	int ndim=0;
	for (int i=0;i<nblk;i++) {
	  comps[i]=new LocalDataContainerSection<DataType>(ldc,ndim,k[i]);
	  ndim+=k[i];
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from PartitionedLocalDataContainer constructor\n";
	throw e;
      }
    }

    /** Main constructor. Creates STL vector of LocalDataContainerSection *'s. */
    PartitionedLocalDataContainer( LocalDataContainer<DataType> & _ldc, 
				   vector<int> _k)
      : ldc(_ldc), k(_k), comps(k.size()) {
      try {
	int nblk = k.size();
	int ndim = 0;
	for (int i=0;i<nblk;i++) {
	  ndim+=k[i];
	}
	if (ndim != ldc.getSize()) {
	  RVLException e;
	  e<<"Error: PartitionedLocalDataContainer constructor\n";
	  e<<"block dimensions do not add up to overall dimensions\n";
	  e<<"overall dimension = "<<ldc.getSize()<<"\n";
	  for (int i=0; i<nblk; i++) {
	    e<<"dimension of block "<<i<<" = "<<k[i]<<"\n";
	  }
	  throw e;
	}
	ndim=0;
	for (int i=0;i<nblk;i++) {
	  comps[i]=new LocalDataContainerSection<DataType>(ldc,ndim,k[i]);
	  ndim+=k[i];
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from PartitionedLocalDataContainer constructor\n";
	throw e;
      }
    }

    ~PartitionedLocalDataContainer() {
      for (int i=0;i<k.size();i++) {
	if (comps[i]) delete comps[i];
      }
    }
  
    /** delegates getSize() */
    int getSize() const { return ldc.getSize(); }

    /** delegates getData() */
    DataType * getData() { return ldc.getData(); }
    DataType const * getData() const { return ldc.getData(); }

    int getNumberOfComponents() const { return k.size(); }
    LocalDataContainer<DataType> & operator[](int i) {
      if (i<0 || i>=(k.size())) {
	RVLException e;
	e<<"Error: PartitionedLocalDataContainer::getLocalComponent\n";
	e<<"input index "<<i<<" out of range [0, "<<(int)(k.size())-1<<"]\n";
	throw e;
      }
      return *(comps[i]);
    }

    LocalDataContainer<DataType> const & operator[](int i) const {
      if (i<0 || i>=(k.size())) {
	RVLException e;
	e<<"Error: PartitionedLocalDataContainer::getLocalComponent\n";
	e<<"input index "<<i<<" out of range [0, "<<(int)(k.size())-1<<"]\n";
	throw e;
      }
      return *(comps[i]);
    }

    ostream & write(ostream & str) const {
      str<<"Partitioned LocalDataContainer Object\n";
      str<<"   number of components = "<<getNumberOfComponents()<<"\n";
      for (int i=0;i<getNumberOfComponents();i++) {
	str<<"   *** component "<<i<<"\n";
	comps[i]->write(str);
      }
      return str;
    }
  };

}

#endif
