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

#ifndef __HCL2RN
#define __HCL2RN

#include "functions.hh"
#include "rn.hh"
#include "localspace.hh"
#include "locallinalg.hh"

namespace RVL {

  /** A space whose elements are RnArrays and which uses the standard
      linear algebra package.  This is an implementation of the mathematical
      space \f$R^{n}\f$.
  */
  template<class Scalar>
  class RnSpace: public LocalSpace<Scalar> {
  private:

    int dim;
    RnDataContainerFactory<Scalar> * dcfac;
    RVLLinearAlgebraPackage<Scalar> lap;
    RVLLinearAlgebraPackage<Scalar> & lapref;

  protected:

    void redim(int newdim) { 
      try {
	dim = newdim;
	if (dcfac) delete dcfac;
	dcfac = new RnDataContainerFactory<Scalar>(dim);
      }
      catch (RVLException & e) {
	e<<"\ncalled from RnSpace::redim\n";
	throw e;
      }
    }

    LocalDataContainerFactory<Scalar> & getLDCF() const { return *dcfac; }
    LinearAlgebraPackage<Scalar> & getLAP() const { return lapref; }

  public:

    RnSpace(int n=0): dim(n), lap(), lapref(lap) {
      try {
	dcfac = new RnDataContainerFactory<Scalar>(dim);
      }
      catch (RVLException & e) {
	e<<"\ncalled from RnSpace constructor\n";
	throw e;
      }
    }
    RnSpace(const RnSpace<Scalar> & sp): dim(sp.dim), lap(), lapref(lap) {
      try {
	dcfac = new RnDataContainerFactory<Scalar>(dim);
      }
      catch (RVLException & e) {
	e<<"\ncalled from RnSpace constructor\n";
	throw e;
      }
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

    virtual ~RnSpace() {
      if (dcfac) delete dcfac;
    }

    int getSize() const { return dim; }

    virtual void write(RVLException & str) const {
      str<<"RnSpace: simple dense vector space\n";
      str<<"dimension = "<<dim<<"\n";
    }
    virtual ostream & write(ostream & str) const {
      str<<"RnSpace: simple dense vector space\n";
      str<<"dimension = "<<dim<<"\n";
      return str;
    }
  };

}

#endif







