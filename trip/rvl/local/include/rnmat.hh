/*************************************************************************

Copyright Rice University, 2011.
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

#ifndef __RVL_RNMAT
#define __RVL_RNMAT

#include "op.hh"
#include "local.hh"
#include "rnspace.hh"

namespace RVL {

  /** Mat-Vec as FunctionObject. Stores matrix as an RnArray. Provides
      action of matrix or of its transpose. */

  template<typename T>
  class matvec:  public BinaryLocalFunctionObject<T> {

  private: 
    
    int rows;
    int cols;
    mutable bool adj;
    RnArray<T> mat;
    matvec();

  public:

    matvec(int _rows, int _cols): rows(_rows), cols(_cols), adj(false), mat(rows*cols) {
      for (int i=0;i<rows*cols; i++) mat.getData()[i]=ScalarFieldTraits<T>::Zero();
    }
    matvec(matvec<T> const * m): rows(m->rows), cols(m->cols), adj(m->adj), mat(m->mat) {}
    ~matvec() {}

    /** expose data */
    T * getData() { return mat.getData(); }
    T const * getData() const { return mat.getData(); }

    /** do something to all elements */
    void eval(UnaryLocalFunctionObject<T> & f) { f(mat); }

    /** mutable element access */
    T & getElement(int i, int j) {
      if (i<0 || i>rows-1) {
	RVLException e;
	e<<"Error: matvec::getElement\n";
	e<<"row index "<<i<<" out of range [0,"<<rows-1<<"]\n";
	throw e;
      }
      if (j<0 || j>cols-1) {
	RVLException e;
	e<<"Error: matvec::getElement\n";
	e<<"col index "<<j<<" out of range [0,"<<cols-1<<"]\n";
	throw e;
      }
      return mat.getData()[i+rows*j];
    }

    /** const element access */
    T const & getElement(int i, int j) const {
      if (i<0 || i>rows-1) {
	RVLException e;
	e<<"Error: matvec::getElement\n";
	e<<"row index "<<i<<" out of range [0,"<<rows-1<<"]\n";
	throw e;
      }
      if (j<0 || j>cols-1) {
	RVLException e;
	e<<"Error: matvec::getElement\n";
	e<<"col index "<<j<<" out of range [0,"<<cols-1<<"]\n";
	throw e;
      }
      return mat.getData()[i+rows*j];
    }

    /** toggle between operator and adjoint */
    void setAdj(bool flag) { adj=flag; }
    /** tell the truth! */
    bool getAdj() const { return adj; }

    /** Note that we use the RVL::conj function defined in utility.hh
	- this delegates to std::complex::conj when appropriate,
	otherwise is no-op, which is correct for all floating point
	data types in std.
    */
    using RVL::BinaryLocalEvaluation<T>::operator();
    void operator()(LocalDataContainer<T> & y, LocalDataContainer<T> const & x) {

      if (adj) {

	if (y.getSize() < (size_t(cols)) || x.getSize() < size_t(rows)) {
	  RVLException e;
	  e<<"Error: matvec::operator(), adjoint\n";
	  e<<"either input or output too short\n";
	  throw e;
	}
	for (int j=0;j<cols;j++) {
	  y.getData()[j]=ScalarFieldTraits<T>::Zero();
	  for (int i=0;i<rows;i++) 
	    y.getData()[j]+= RVL::conj(mat.getData()[i+j*rows])*x.getData()[i];
	}
	
      }
      else {

	if (y.getSize() < size_t(rows) || x.getSize() < size_t(cols)) {
	  RVLException e;
	  e<<"Error: matvec::operator()\n";
	  e<<"either input or output too short\n";
	  e<<"input size = "<<x.getSize()<<" should be at least "<<cols<<"\n";
	  e<<"output size = "<<y.getSize()<<" should be at least "<<rows<<"\n";
	  throw e;
	}
	for (int i=0;i<rows;i++) {
	  y.getData()[i]=ScalarFieldTraits<T>::Zero(); 
	  for (int j=0;j<cols;j++) 
	    y.getData()[i]+=mat.getData()[i+j*rows]*x.getData()[j];
	}
      }
    }

    string getName() const { string tmp = "matvec"; return tmp; }

  };

  /** forward action wrapper - for generation of LinearOpFO instances */
  template<typename T>
  class fmatvec:  public BinaryLocalFunctionObject<T> {
    
  private:
    matvec<T> & m;
    fmatvec();

  public:
    fmatvec( matvec<T> & _m ) : m(_m) {}
    fmatvec( fmatvec<T> const & f) : m(f.m) {}
    ~fmatvec() {}

    using RVL::BinaryLocalEvaluation<T>::operator();
    void operator()(LocalDataContainer<T> & y,
		    LocalDataContainer<T> const & x) {
      m.setAdj(false);
      m(y,x);
    }
    string getName() const { string tmp = "fmatvec"; return tmp; }
  };

  /** adjoint action wrapper - for generation of LinearOpFO instances */
  template<typename T>
  class amatvec:  public BinaryLocalFunctionObject<T> {
    
  private:
    matvec<T> & m;
    amatvec();

  public:
    amatvec( matvec<T> & _m ) : m(_m) {}
    amatvec( amatvec<T> const & f) : m(f.m) {}
    ~amatvec() {}

    using RVL::BinaryLocalEvaluation<T>::operator();
    void operator()(LocalDataContainer<T> & y,
		    LocalDataContainer<T> const & x) {
      m.setAdj(true);
      m(y,x);
    }
    string getName() const { string tmp = "amatvec"; return tmp; }
  };

  /** GenMat: simple matrix class, mostly for testing purposes. A
      serious Matrix library would build on one of the standard
      libraries, like GSL or even LAPACK. This set of classes is
      self-contained and sufficient for testing and simple
      low-dimensional apps. */

  template<typename T>
  class GenMat: public LinearOp<T> {

  private: 

    RnSpace<T> dom;
    RnSpace<T> rng;
    mutable matvec<T> a;
    GenMat();

  protected: 

    virtual LinearOp<T> * clone() const { return new GenMat(*this); }

    void apply(Vector<T> const & x,
	       Vector<T> & y) const {
      try {
	a.setAdj(false);
	y.eval(a,x);
      }
      catch (RVLException & e) {
	e<<"\ncalled from GenMat::apply\n";
	throw e;
      }
    }

    void applyAdj(Vector<T> const & x,
		  Vector<T> & y) const {
      try {
	a.setAdj(true);
	y.eval(a,x);
      }
      catch (RVLException & e) {
	e<<"\ncalled from GenMat::applyAdj\n";
	throw e;
      }
    }

  public:

    /** main constructor requires only domain and range */
    GenMat(RnSpace<T> const & _dom, RnSpace<T> const & _rng)
      : dom(_dom), rng(_rng), a(rng.getSize(),dom.getSize()) {}

    /** copy constructor, used in clone */
    GenMat(GenMat<T> const & m) 
      : dom(m.dom), rng(m.rng), a(m.a)  {}

    ~GenMat() {}

    Space<T> const & getDomain() const { return dom; }
    Space<T> const & getRange() const { return rng; }

    /** convenience access to row and col dims */
    int getNRows() const { return int(rng.getSize()); }
    int getNCols() const { return int(dom.getSize()); }

    /** access to data - essential for typical apps, eg. calls to LAPACK */
    T * getData() { return a.getData(); }
    T const * getData() const { return a.getData(); }

    /** element access - not so useful, but bound-checked and why not */
    // const
    T const & getElement(int i, int j) const { 
      try { return a.getElement(i,j); }
      catch (RVLException e) {e<<"\ncalled from GenMat::getElement\n"; throw e; }
    }

    // mutable
    T & getElement(int i, int j) {
      try { return a.getElement(i,j); }
      catch (RVLException e) {e<<"\ncalled from GenMat::getElement\n"; throw e; }
    }

    /* set - could be used for bound-checked assignemnt */
    virtual void setElement(int i, int j, T e) {
      try { a.getElement(i,j)=e; }
      catch (RVLException e) {e<<"\ncalled from GenMat::setElement\n"; throw e; }
    }

    /* do something to the the internal LDC via a FO eval */
    void eval(UnaryLocalFunctionObject<T> & f) { a.eval(f); }

    virtual ostream & write(ostream & str) const {
      str<<"GenMat: simple general matrix class\n";
      str<<"based on matvec FunctionObject\n";
      str<<"rows = "<<rng.getSize()<<" cols = "<<dom.getSize()<<endl;
      return str;
    }
    
  };
  
  /** Symmetric specialization of GenMat. Differs only in element
      assignment, which forces symmetry */
  template<typename T>
  class SymMat: public GenMat<T> {
    
  protected:
    
    LinearOp<T> * clone() const { return new SymMat<T>(*this); }

  public:

    SymMat(RnSpace<T> const & dom)
      : GenMat<T>(dom,dom) {}

    SymMat(SymMat<T> const & m) 
      : GenMat<T>(m) {}

    ~SymMat() {}

    void setElement(int i, int j, T e) {
      try { GenMat<T>::setElement(i,j,e); GenMat<T>::getElement(j,i,e); }
      catch (RVLException e) {e<<"\ncalled from GenMat::setElement\n"; throw e; }
    }

    ostream & write(ostream & str) const {
      str<<"SymMat: simple symmetric matrix class\n";
      str<<"based on matvec FunctionObject\n";
      RnSpace<T> const & dom = dynamic_cast<RnSpace<T > const &>(this->getDomain());
      str<<"rows = cols = "<<dom.getSize()<<endl;
      return str;
    }
  };

}

#endif
