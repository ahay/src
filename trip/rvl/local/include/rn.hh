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

#ifndef __RVLRN
#define __RVLRN

#include "local.hh"

namespace RVL {

  /** An implementation of a LDC containing an array of Scalars.
      This class will likely serve as the base class for
      a wide variety of LDCs.
  */
  template<class Scalar>
  class RnArray: public LocalDataContainer<Scalar> {
  private:

    size_t n;
    int own;
    size_t start;
    mutable Scalar * a;

    RnArray() {}

  public:

    /** deep copy, following example of stl containers */
    RnArray(const RnArray<Scalar> & x)
      : n(x.n), own(1), start(0), a(new Scalar[n]) {
      for (size_t i=0;i<n;i++) {
	a[i]=x.a[i];
      }
    }

    /** main constructor: dimension passed as sole argument. */
    RnArray(size_t _n)
      : n(_n),
	own(1), start(0), a(NULL) {
      if (n<=0) {
	RVLException e;
	e<<"Error: RnArray constructor\n";
	e<<"must have positive dimension = "<<n<<"\n";
	throw e;
      }
    }

    /** constructor which gives a "view" of another RnArray. */
    RnArray(LocalDataContainer<Scalar> & rn, size_t len, size_t _start=0)
      : start(_start), n(len), own(0) {
      a = &(rn.getData()[start]);
      if (start+n>rn.getSize()) {
	RVLException e;
	e<<"Error: RnArray \"view\" constructor\n";
	e<<"segment runs off end of input LocalDC\n";
	e<<"length of input LocalDC = "<<rn.getSize()<<"\n";
	e<<"index of segment start  = "<<start<<"\n";
	e<<"length of segment       = "<<n<<"\n";
	throw e;
      }
    }

    ~RnArray() { if (own && a) delete [] a; }

    size_t getSize() const { return n; }

    Scalar * getData() { 
      if (!a) a = new Scalar[n];
      return a; 
    }

    Scalar const * getData() const {
      if (!a) a = new Scalar[n];
      return a;
    }

    void write(RVLException & str) const {
      str<<"RnArray Local Data Container object\n";
      str<<"length = "<<n<<"\n";
      if (a) {
	str<<"samples: \n";
	for (size_t i=0;i<n;i++) {
	  str<<"data["<<i<<"] = "<<a[i]<<"\n";
	}
      }
      else {
	str<<"data samples not initialized\n";
      }
    }

    ostream & write(ostream & str) const {
      str<<"RnArray Local Data Container object\n";
      str<<"length = "<<n<<"\n";
      str<<"samples: \n";
      if (a) {
	for (size_t i=0;i<n;i++) {
	  str<<"data["<<i<<"] = "<<a[i]<<"\n";
	}
      }
      else {
	str<<"data samples not initialized\n";
      }
      return str;
    }
  };

  /** A factory for RnArrays.  */
  template<class Scalar> 
  class RnDataContainerFactory: public LocalDataContainerFactory<Scalar> {

  private:

    size_t n;

  protected:

    void redim(size_t _n) { n = _n; }

  public:

    RnDataContainerFactory(size_t _n=0): n(_n) {}
    RnDataContainerFactory(const RnDataContainerFactory<Scalar> & f): n(f.n) {}
    ~RnDataContainerFactory() {}

    virtual LocalDataContainer<Scalar> * buildLocal() const { 
      try {
	return new RnArray<Scalar>(n);
      }
      catch (RVLException & e) {
	e<<"\ncalled from RnDataContainerFactory::buildLocal()\n";
	throw e;
      }
    }      

    /** determines whether another DataContainerFactory builds the same
	kind of DataContainer. Usually this means effectively that the two
	are of the same type, determined by runtime type checking. Returns
	zero if different types, otherwise nonzero.
    */
    virtual bool compare( DataContainerFactory const & dcf) const {
      try {
	if (&dcf == this) return 1;
	RnDataContainerFactory<Scalar> const & a = 
	  dynamic_cast< RnDataContainerFactory<Scalar> const &>(dcf);
	if (n == a.n) return 1;
	return 0;
      }
      catch (bad_cast) {
	return 0;
      }
    }

    /** determines whether a DataContainer is of the type built by this
	factory. Usually implemented through runtime type checking.
	Returns zero if not, else nonzero. 
    */
    virtual bool isCompatible(DataContainer const & dc) const  {
      try {
	RnArray<Scalar> const & a = 
	  dynamic_cast<RnArray<Scalar> const &>(dc);
	if (n == a.getSize()) return 1;
	return 0;
      }
      catch (bad_cast) {
	return 0;
      }
    }

    size_t getSize() const { return n; }

    virtual void write(RVLException & e) const {
      e<<"RnDataContainerFactory, dimn = "<<n<<"\n";
    }
    virtual ostream & write(ostream & str) const {
      str<<"RnDataContainerFactory, dimn = "<<n<<"\n";
      return str;
    }

  };

}

#endif







