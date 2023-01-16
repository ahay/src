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

#ifndef __RVL_PVEC
#define __RVL_PVEC

#include "space.hh"
#include "productdata.hh"

namespace RVL {

  /** Abstract base class expressing behaviour of Cartesian products of
      vector spaces. Uses the ProductDataContainer class to define data
      component access. Provides basic data container methods (getSize,
      operator[]) for a Cartesian product. All other methods of Space
      are implemented in terms of these. 

      A default implementation is provided for inner product, expressing
      the canonical inner product in a product Hilbert space. This method
      is virtual and may be overridden in subclasses.  */

  template<class Scalar>
  class ProductSpace: public Space<Scalar>,
		      public ROProduct< Space<Scalar> > {

  public:

    ProductSpace() {}
    ProductSpace(const ProductSpace<Scalar> &) {}
#if __cplusplus >= 201103L // We are using C++11 or a later version
    virtual ~ProductSpace() noexcept(false) {}
#else
    virtual ~ProductSpace() {}
#endif
    virtual bool isCompatible(DataContainer const & dc) const {
      try {
	ProductDataContainer const & pdc =
	  dynamic_cast<ProductDataContainer const &>(dc);
	if (pdc.getSize() != this->getSize()) return 0;
	for (size_t i=0;i<this->getSize();i++) {
	  if (!((*this)[i].isCompatible(pdc[i]))) return 0;
	}
	return true;
      }
      catch (bad_cast const&) {
	return false;
      }
    }

    virtual bool operator ==(const Space<Scalar> & sp) const {
      if (this == &sp) return true;
      try {
	const ProductSpace<Scalar> & psp 
	  = dynamic_cast<const ProductSpace<Scalar> &>(sp);
	if (psp.getSize() != this->getSize()) return false;
	for (size_t i=0;i<this->getSize();i++) {
	  if (psp[i] != (*this)[i]) return false;
	}
	return true;
      }
      catch (bad_cast const&) {
	return false;
      }
    }

    /** inner product - may be overridden by supplying non-block-diag
	Gram matrix. */
    virtual Scalar inner(DataContainer const & x, 
			 DataContainer const & y) const {
      Scalar res = ScalarFieldTraits<Scalar>::Zero();
      try {
	ProductDataContainer const & xp = 
	  dynamic_cast<ProductDataContainer const &>(x);
	ProductDataContainer const & yp = 
	  dynamic_cast<ProductDataContainer const &>(y);
	if (xp.getSize() != this->getSize() || yp.getSize() != this->getSize()) {
	  RVLException e; e<<"Error: ProductSpace::inner\n";
	  e<<"input data containers not of same size as space\n";
	  throw e;
	}
	for (size_t i=0;i<this->getSize();i++) {
	  try {
	    res += (*this)[i].inner(xp[i],yp[i]);
	  }
	  catch (RVLException & e) {
	    e<<"\ncalled from ProductSpace::inner\n";
	    throw e;
	  }
	}
	return res;
      }
      catch (bad_cast const&) {
	RVLException e; e<<"Error: ProductSpace::inner\n";
	e<<"inputs not (both) product vectors\n";
	throw e;
      }
    }

    // zero vector
    virtual void zero(DataContainer & x) const {
      try {
	ProductDataContainer & xp = 
	  dynamic_cast<ProductDataContainer &>(x);
	if (xp.getSize() != this->getSize()) {
	  RVLException e; e<<"Error: ProductSpace::zero\n";
	  e<<"input data container not of same size as space\n";
	  throw e;
	}
	for (size_t i=0;i<this->getSize();i++) {
	  (*this)[i].zero(xp[i]);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from ProductSpace::zero\n";
	throw e;
      }
      catch (bad_cast const&) {
	RVLException e; e<<"Error: ProductSpace::zero\n";
	e<<"input not product vector\n";
	throw e;
      }
    }

    // linear combination
    virtual void linComb(Scalar a, DataContainer const & x,
			 Scalar b, DataContainer & y) const {
      try {
	ProductDataContainer const & xp = 
	  dynamic_cast<ProductDataContainer const &>(x);
	ProductDataContainer & yp = 
	  dynamic_cast<ProductDataContainer &>(y);
	if (xp.getSize() != this->getSize() || 
	    yp.getSize() != this->getSize()) {
	  RVLException e; e<<"Error: ProductSpace::linComb\n";
	  e<<"input data containers not of same size as space\n";
	  e<<"**** FIRST INPUT PDC:\n";
	  x.write(e);
	  e<<"**** SECOND INPUT PDC:\n";
	  y.write(e);
	  e<<"**** SPACE:\n";
	  Space<Scalar>::write(e);
	  throw e;
	}
	for (size_t i=0;i<this->getSize();i++) {
	  (*this)[i].linComb(a,xp[i],b,yp[i]);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from ProductSpace::linComb\n";
	throw e;
      }
      catch (bad_cast const&) {
	RVLException e; e<<"Error: ProductSpace::linComb\n";
	e<<"inputs not all product vectors\n";
	throw e;
      }
    }

    /** Write method - may be overwritten if useful. */
    virtual ostream & write(ostream & str) const {
      str<<"Product Space"<<endl;
      str<<"size = "<<this->getSize()<<"\n";
      for (size_t i=0; i<this->getSize(); i++) {
	str<<"***factor "<<i<<":"<<endl;
	(*this)[i].write(str);
      }
      return str;
    }
  };

  /** Standard construction of product spaces, via an (STL) vector of
      Spaces. Note that all data references are to external objects. 
      Therefore copy is shallow. No data other than pointers is owned.

      In order avoid dangling pointer issues, the only usable
      constructors take explicit lists of const space references. Thus
      StdProductSpace can be instantiated only within the scope of the
      component Space instantiations.

      Convenience constructors supplied for 2, 3, and 4
      components. Arbitrary products may be built up via the
      associativity of Cartesian product.
  */

  template<class Scalar>
  class StdProductSpace: public ProductSpace<Scalar> {

  private: 

    /** vector of const pointers to Space */
    std::vector< Space<Scalar> const * > s;

    /** Default construction - disallowed */
    StdProductSpace();

  public:

    /** Copy construction. [should this be public? WWS 16.10.12 */
    StdProductSpace(const StdProductSpace<Scalar> & sp): s(sp.s.size()) {
      for (size_t i=0;i<s.size();i++) {
	s[i]=(sp.s)[i];
      }
    }

    /** constructor for single space - allows one to regard a single space
	as a (trivial) product space */
    StdProductSpace(Space<Scalar> const & s1): s() {
      s.push_back(&s1);
    }

    /** constructor for pairs of spaces. */
    StdProductSpace(Space<Scalar> const & s1,
		    Space<Scalar> const & s2): s() {
      s.push_back(&s1);
      s.push_back(&s2);
    }
    /** constructor for triples of spaces. */
    StdProductSpace(Space<Scalar> const & s1,
		    Space<Scalar> const & s2,
		    Space<Scalar> const & s3): s() {
      s.push_back(&s1);
      s.push_back(&s2);
      s.push_back(&s3);
    }
    /** constructor for quadruples of spaces. */
    StdProductSpace(Space<Scalar> const & s1,
		    Space<Scalar> const & s2,
		    Space<Scalar> const & s3,
		    Space<Scalar> const & s4): s() {
      s.push_back(&s1);
      s.push_back(&s2);
      s.push_back(&s3);
      s.push_back(&s4);
    }

    /** Destructor does not destroy anything except the STL vector
	structure (and that is done automatically), as all references
	are external. */
    ~StdProductSpace() {}

    /** implements virtual DataContainer constructor via
	StdProductDataContainer class. */
    DataContainer * buildDataContainer() const {
      StdProductDataContainer * d = 
	new StdProductDataContainer();
      for (size_t i=0; i<getSize(); i++) {
	SpaceDCF<Scalar> f(*(s[i]));
	d->push(f);
      }
      return d;
    }

    size_t getSize() const { return s.size(); }
    Space<Scalar> const & operator[](size_t i) const { return *(s[i]); }

  };


  /**  This space implements the Cartesian power of a space.
       It requires the base space and the power as input,
       then implements all the Space methods.
       
       This is a child of ProductSpace, since it makes 
       logical sense to be so, and it can be compared
       to other ProductSpaces.  However, the methods
       have been reimplemented from ProductSpace
       in order to avoid inefficiency.
  */
  template<class Scalar>
  class CartesianPowerSpace: public ProductSpace<Scalar> {

  protected:
    
    size_t size;
    const Space<Scalar> & subspc;

    CartesianPowerSpace() {}
    
  public:
    CartesianPowerSpace(size_t _size, const Space<Scalar> &  _subspc)
      :size(_size), subspc(_subspc) {}
    
    CartesianPowerSpace(const CartesianPowerSpace<Scalar> & s)
      :size(s.size), subspc(s.subspc) {}
    
    virtual ~CartesianPowerSpace() {}
    
    /** returns number of components */
    size_t getSize() const { return size; }
    
    /** returns ith component space */
    Space<Scalar> const & operator[](size_t i) const {
      return subspc;
    }

    bool isCompatible(DataContainer const & dc) const {
      try {
	ProductDataContainer const & pdc =
	  dynamic_cast<ProductDataContainer const &>(dc);
	if (pdc.getSize() != getSize()) return 0;
	for (size_t i=0;i<getSize();i++) {
	  if (!(subspc.isCompatible(pdc[i]))) return 0;
	}
	return true;
      }
      catch (bad_cast const&) {
	return false;
      }
    }
    
    bool operator ==(const Space<Scalar> & sp) const {
      if (this == &sp) return true;
      try {
	const CartesianPowerSpace<Scalar> & psp 
	  = dynamic_cast<const CartesianPowerSpace<Scalar> &>(sp);
	if (psp.getSize() != getSize()) return false;
	if (psp[0] != subspc) return false;
	return true;
      }
      catch (bad_cast const&) {
	try {
	  const ProductSpace<Scalar> & psp 
	    = dynamic_cast<const ProductSpace<Scalar> &>(sp);	
	  if (psp.getSize() != getSize()) return 0;
	  for (size_t i=0;i<getSize();i++) {
	    if (subspc != psp[i]) return 0;
	  }
	  return true;
	}
	catch (bad_cast const&) {
	  return false;
	}
      }
    }
    
    // inner product
    Scalar inner(DataContainer const & x, 
		 DataContainer const & y) const {
      Scalar res = 0.0;
      try {
	ProductDataContainer const & xp = 
	  dynamic_cast<ProductDataContainer const &>(x);
	ProductDataContainer const & yp = 
	  dynamic_cast<ProductDataContainer const &>(y);
	if (xp.getSize() != getSize() || yp.getSize() != getSize()) {
	  RVLException e; e<<"Error: CartesianPowerSpace::inner\n";
	  e<<"input data containers not of same size as space\n";
	  throw e;
	}
	for (size_t i=0;i<getSize();i++) {
	  try {
	    res += subspc.inner(xp[i],yp[i]);
	  }
	  catch (RVLException & e) {
	    e<<"\ncalled from CartesianPowerSpace::inner\n";
	    throw e;
	  }
	}
	return res;
      }
      catch (bad_cast const&) {
	RVLException e; e<<"Error: CartesianPowerSpace::inner\n";
	e<<"inputs not (both) product vectors\n";
	throw e;
      }
    }
    
    // zero vector
    void zero(DataContainer & x) const {
      try {
	ProductDataContainer & xp = 
	  dynamic_cast<ProductDataContainer &>(x);
	if (xp.getSize() != getSize()) {
	  RVLException e; e<<"Error: CartesianPowerSpace::zero\n";
	  e<<"input data container not of same size as space\n";
	  throw e;
	}
	for (size_t i=0;i<getSize();i++) {
	  subspc.zero(xp[i]);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from CartesianPowerSpace::zero\n";
	throw e;
      }
      catch (bad_cast const&) {
	RVLException e; e<<"Error: CartesianPowerSpace::zero\n";
	e<<"input not product vector\n";
	throw e;
      }
    }
    
    // linear combination
    void linComb(Scalar a, DataContainer const & x,
		 Scalar b, DataContainer & y) const {
      try {
	ProductDataContainer const & xp = 
	  dynamic_cast<ProductDataContainer const &>(x);
	ProductDataContainer & yp = 
	  dynamic_cast<ProductDataContainer &>(y);
	if (xp.getSize() != getSize() || 
	    yp.getSize() != getSize()) {
	  RVLException e; e<<"Error: CartesianPowerSpace::linComb\n";
	  e<<"input data containers not of same size as space\n";
	  throw e;
	}
	for (size_t i=0;i<getSize();i++) {
	  subspc.linComb(a,xp[i],b,yp[i]);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from CartesianPowerSpace::linComb\n";
	throw e;
      }
      catch (bad_cast const&) {
	RVLException e; e<<"Error: CartesianPowerSpace::linComb\n";
	e<<"inputs not all product vectors\n";
	throw e;
      }
    }
    
    /** implements virtual DataContainer constructor via
	StdProductDataContainer class. */
    DataContainer * buildDataContainer() const {
      StdProductDataContainer * d = 
	new StdProductDataContainer();
      SpaceDCF<Scalar> f(subspc);
      for (size_t i=0; i<getSize(); i++) d->push(f);
      return d;
    }
    
    ostream & write(ostream & str) const {
      str<<"CartesianPowerSpace"<<endl;
      str<<"size = "<<getSize()<<"\n";
      str <<"based on subspace : "<< endl;
      subspc.write(str);
      return str;
    }
    
  };
  
  /** Provides indexed access to components (factors) of a Vector in 
      ProductSpace.  
      
      Set up to work properly for non-product vectors also - these have
      one component!
  */
  template<class Scalar>
  class Components: public Product< Vector<Scalar> > {

  private:

    size_t size;
    int newflag;

    ProductDataContainer * pdc;
    const ProductSpace<Scalar> * psp;

    vector< Vector<Scalar>* > comp;

    Components();
    Components(const Components<Scalar> &);

  public:

    /* Constructor. Uses protected Vector constructor to build
       (STL) vector of component Vector(s). */
    Components(const Vector<Scalar> & v)
	: size(1), pdc(NULL), psp(NULL), comp(1) {
      if ((pdc = 
	   dynamic_cast<ProductDataContainer *>
	   (v.getDataContainer())) 
	  &&
	  (psp = dynamic_cast<const ProductSpace<Scalar> *>(&(v.getSpace())))) {
	if (pdc->getSize() == psp->getSize()) {
	  size = psp->getSize();
	  comp.resize(size);
	  for( size_t i=0;i<size;i++ ) {
	    comp[i]=
	      new Vector<Scalar>((*psp)[i],&((*pdc)[i]),
				 v.getVersionRef(), false);
	  }
	  newflag=1;
	}
      }
      else {
	comp[0]=const_cast<Vector<Scalar> *>(&v);
	newflag=0;
      }
    }

    /** Destructor. Deletes any owned objects. */
    ~Components() {
      if( newflag )
	for( size_t i=0;i<size;i++ ) delete comp[i];
    }

    /** return number of components */
    size_t getSize() const { return size; }
    /** return ith component; throw exception if index is out of bounds. */
    Vector<Scalar> & operator[](size_t i) {
      if (i>size-1) {
	RVLException e; e<<"Error: Components::operator[]\n"; 
	e<<"index "<<i<<" out of range [0, "<<size-1<<"]\n";
	throw e;
      }
      return *(comp[i]);
    }
    Vector<Scalar> const & operator[](size_t i) const {
      if (i>size-1) {
	RVLException e; e<<"Error: Components::operator[]\n"; 
	e<<"index out of range [0, "<<size-1<<"]\n";
	throw e;
      }
      return *(comp[i]);
    }

    ostream & write(ostream & str) const {
      str<<"Components of a Vector"<<endl;
      str<<"  number = "<<size<<endl;
      for (size_t i=0;i<size;i++) {
	str<<"***component "<<i<<":"<<endl;
	(*this)[i].write(str);
      }
      return str;
    }
  };
}

#endif
