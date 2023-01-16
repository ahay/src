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

#ifndef __RVL_PDC
#define __RVL_PDC

#include "data.hh"
#include "product.hh"

namespace RVL {

  /** BlockFOs are arrays of FOs, which can be evaluated
      block-diagonal fashion on ProductDCs. Production of BlockFORs is
      more involved, put off until need arises. Concrete implementation 
      using std::vector.
  */
  class BlockFunctionObject: 
    public FunctionObject,
    public std::vector<FunctionObject *> {

  public:
    
    string getName() const {
      string tmp;
      tmp+="*****************************\n";
      tmp+="*** Block Function Object ***\n";
      tmp+="*****************************\n";
      /*
      tmp+="size = ";
      tmp+=this->size()<<"\n";
      for (int i=0;i<this->size();i++) 
	tmp<<"component "<<i<<" name = "<<this->at(i).getName()<<"\n";
	tmp<<"*****************************\n";
      */
      return tmp;
    }
  };

  /** Produces diagonal ("scalar") app of single FO to all components
      of ProductDC. 
  */
  class DiagonalFunctionObject: public BlockFunctionObject {
  private:
    DiagonalFunctionObject();
    DiagonalFunctionObject(DiagonalFunctionObject const &);

  public:
    DiagonalFunctionObject(size_t n, FunctionObject & f)
      : BlockFunctionObject() {
      for (size_t i=0;i<n;i++) this->push_back(&f);
    }
  };

  /** ProductDataContainers are DataContainers equipped with an
      indexing operator[], which returns a reference to a
      DataContainer when supplied with an in-range int index.
      Since ProductDataContainers act in this way like arrays of
      DataContainers, they also have a getSize() method. Both the
      indexing operator and getSize() are supplied by a mixin Product
      interface.

      This class is abstract, to permit a variety of schemes for
      storing and retrieving references to the component
      DataContainers. The non-const evaluation methods inherited from
      DataContainer are designed for BlockFOs, that is, effectively
      matrices of FOs. An option is provided to create a DiagonalFO;
      in that case evaluation is effectively equivalent to the simple
      loop algorithms of the Standard Library, but defined at a more
      ablstract level, not requiring copy semantics of the items over
      which the iteration takes place. For FOCEs, only the diagonal
      option is provided: evaluation proceeds by sequential evaluation
      of the FOCE on the components. To override this assumption, it
      will be necessary to view the object as another type of
      DataContainer (i.e. other than ProductDataContainer) and invoke
      an alternative implementation of evaluation.
  */

  class ProductDataContainer: public DataContainer, 
			      public Product<DataContainer> {

  public:

    ProductDataContainer() {}
    ProductDataContainer(const ProductDataContainer &) {}
#if __cplusplus >= 201103L  // We are using C++11 or a later version
    virtual ~ProductDataContainer() noexcept(false) {}
#else
    virtual ~ProductDataContainer() {}
#endif
    void eval( FunctionObject & f,
	       std::vector<DataContainer const *> & x) {
      /* first create block FO either by cast or diag construction */
      BlockFunctionObject * bf = NULL;
      if (!(bf=dynamic_cast<BlockFunctionObject *>(&f))) {
	bf=new DiagonalFunctionObject(this->getSize(),f);
      }
      if (!bf) {
	RVLException e;
	e<<"Error: ProductDataContainer::eval\n";
	e<<"failed to create BlockFO from FO = "<<f.getName()<<" by either cast or DiagFO\n";
	throw e;
      }	
      if (bf->size()!=this->getSize()) {
	RVLException e;
	e<<"Error: ProductDataContainer::eval\n";
	e<<"Input FO = "<<f.getName()<<" cast to BlockFO but wrong size\n";
	e<<"FO size = "<<bf->size()<<" ProductDC size = "<<this->getSize()<<"\n";
	throw e;
      }

      try {
	size_t nx = x.size();
	vector<ProductDataContainer const *> xp(nx);
	for (size_t i=0;i<nx;i++) {
	  if (!(xp[i]=dynamic_cast<ProductDataContainer const *> (x[i]))) {
	    RVLException e;
	    e<<"Error: ProductDataContainer::eval\n";
	    e<<"argument "<<i<<" is not PDC\n";
	    throw e;
	  }
	  if (xp[i]->getSize() != getSize()) {
	    RVLException e;
	    e<<"Error: ProductDataContainer::eval \n";
	    e<<"input ProductDataContainer arg "<<i
	     <<" does not have same number\n";
	    e<<"of factors\n";
	    throw e;
	  }      
	}
	vector<DataContainer const *> tmp(nx);
	size_t nc = this->getSize();
	for (size_t i=0;i<nc;i++) {
	  for (size_t j=0;j<nx;j++) {
	    tmp[j] = &((*xp[j])[i]);
	  }
	  //	  (*this)[i].eval(f,tmp);
	  (*this)[i].eval(*(bf->at(i)),tmp);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from ProductDataContainer::eval\n";
	throw e;
      }
      catch (...) {
	throw;
      }
    }

    void eval( FunctionObjectConstEval & f,
	       std::vector<DataContainer const *> & x) const {
      try {
	size_t nx = x.size();
	vector<ProductDataContainer const *> xp(nx);
	for (size_t i=0;i<nx;i++) {
	  if (!(xp[i]=dynamic_cast<ProductDataContainer const *> (x[i]))) {
	    RVLException e;
	    e<<"Error: ProductDataContainer::eval\n";
	    e<<"argument "<<i<<" is not PDC\n";
	    throw e;
	  }
	  if (xp[i]->getSize() != getSize()) {
	    RVLException e;
	    e<<"Error: ProductDataContainer::eval \n";
	    e<<"input ProductDataContainer arg "<<i
	     <<" does not have same number\n";
	    e<<"of factors\n";
	    throw e;
	  }      
	}
	vector<DataContainer const *> tmp(nx);
	size_t nc = this->getSize();
	for (size_t i=0;i<nc;i++) {
	  for (size_t j=0;j<nx;j++) {
	    tmp[j] = &((*xp[j])[i]);
	  }
	  DataContainer const & thisref = (*this)[i];
	  thisref.eval(f,tmp);
	}
      }
      catch (RVLException & e) {
	e<<"\ncalled from ProductDataContainer::eval (const)\n";
	throw e;
      }
      catch (...) {
	throw;
      }
    }

    /** report to ostream */
    ostream & write(ostream & str) const {
      str<<"Product Data Container"<<endl;
      for (size_t i=0; i<getSize(); i++) {
	str<<"***factor "<<i<<":"<<endl;
	(*this)[i].write(str);
      }
      return str;
    }

  };

  /** Standard implementation of ProductDataContainer. Stores
      components by pointer in a std::vector. The default constructor
      creates a StdProductDataContainer which does not own any data,
      and pointers to DataContainers can be appended to it using the
      tt push() method.  

      Major change, WWS, 09.09: no longer a "stupid pointer" class,
      but an intrusive handle like other RVL containers. Instead of a
      pointer to DC, push method accepts a const DCF reference and
      uses the DCF's build method to add a pointer. Thus all data
      referenced by a StdPDC object is managed by it.
  */

  class StdProductDataContainer: 
    public ProductDataContainer,
    public vector< DataContainer * > {

  private:

    /** Copy constructor. Deep copy - creates independent component
	objects, but does not copy data. Really a clone method. 
	Implementation removed due to dependence on DataContainer::clone(),
	which has been deprecated
    */
    StdProductDataContainer(const StdProductDataContainer & p);

  public:

    /** Default constructor yields size=0 object */
    //    StdProductDataContainer(): d() {}
    StdProductDataContainer() {}

    ~StdProductDataContainer() {
      //      for (int i=0;i<getSize();i++) delete d[i];
      for (size_t i=0;i<getSize();i++)
	if (this->at(i)) delete this->at(i);
    }

    size_t getSize() const { return this->size(); }
  
    DataContainer & operator[](size_t i) {
      try {
	return *(this->at(i));
      }
      catch (out_of_range const&) {
	RVLException e;
	e<<"attempt to access component "<<i<<" of StdProductDC of size "<<this->size()<<"\n";
	throw e;
      }
    }
  
    DataContainer const & operator[](size_t i) const {
      try {
	return *(this->at(i));
      }
      catch (out_of_range const&) {
	RVLException e;
	e<<"attempt to access component "<<i<<" of StdProductDC of size "<<this->size()<<"\n";
	throw e;
      }
    }  

    void push(DataContainerFactory const & f) {
      DataContainer * x = f.build();
      this->push_back(x);
    }
  
  };

}

#endif










