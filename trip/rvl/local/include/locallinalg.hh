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

#ifndef __RVL_LLA
#define __RVL_LLA

#include "local.hh"
#include "linalg.hh"
#include "functions.hh"

namespace RVL {

  /** LocalLinCombObject implements linear combination for
      LocalDataContainers: a call to operator() with
      LocalDataContainer arguments y and x in that order should
      implement y <- a*x + b*y, with whatever error checking is
      appropriate. Scalars a and b specified through a setScalar
      method, the only additional interface provided. Ordering is: a
      multiplies const input arg, b multiplies mutable input/output
      arg on rhs.
  */

  /** Interface to local function objects defining the basic ops of
      linear algebra in Hilbert space: linear combination, assignment
      to the zero vector, and inner product. Virtual base implementing
      the LinearAlgebraPackage interface at the local level,
      permitting any convenient definition of these objects.
  */

  template<class DataType, class Scalar = DataType>
  class LocalLinearAlgebraPackage: public LinearAlgebraPackage<Scalar> {

  public:
    LocalLinearAlgebraPackage() {}
    LocalLinearAlgebraPackage(const LocalLinearAlgebraPackage<Scalar, DataType> &) {}
    virtual ~LocalLinearAlgebraPackage() {}
  
    /** access to inner product */
    virtual BinaryLocalFunctionObjectScalarRedn<DataType,Scalar> & localinner() const = 0;
    /** inherited access to inner product */
    FunctionObjectScalarRedn<Scalar> & inner() const { return localinner(); }
    /** access to zero assignment */
    virtual UnaryLocalFunctionObject<DataType> & localzero() const = 0;
    /** inherited access to zero assignment */
    FunctionObject & zero() const { return localzero(); }
    // lincombobject deferred - must be specialized FO type

    /** Compare for compatibility with another LinearAlgebraPackage.
	Usual comparison basis - is the type the same? However
	"compatibility" can be defined more loosely when appropriate:
	the intended meaning is "produces the same results when
	applied to the same data". Returns zero if not compatible,
	nonzero otherwise. 
    */
    virtual bool compare( LinearAlgebraPackage<Scalar> const & lap) const = 0;

    /** report to exception */
    virtual void write(RVLException & e) const = 0;

    /** report to ostream */
    virtual ostream & write(ostream & str) const = 0;
  };

  /** This binary function object performs a linear combination
      u = a*v + b*u. Standard implementation, should be usable across
      a wide variety of data structures underlying vector algebra.     
  */
  template<class Scalar>
  class RVLLinCombObject: public BinaryLocalEvaluation<Scalar>, public LinCombObject<Scalar> {
  private:
    mutable Scalar a;
    mutable Scalar b;
  public:
    RVLLinCombObject(Scalar ain = ScalarFieldTraits<Scalar>::One(), 
		     Scalar bin = ScalarFieldTraits<Scalar>::One()): a(ain), b(bin) {}
    RVLLinCombObject(const RVLLinCombObject<Scalar> & lc ): a(lc.a), b(lc.b) {}
    ~RVLLinCombObject() {}

    /** runtime assignment of coefficients */
    void setScalar(Scalar ain, Scalar bin) { a=ain; b=bin; }

    /** u = a*v + b*u */
    using RVL::BinaryLocalEvaluation<Scalar>::operator();
    void operator() (LocalDataContainer<Scalar> & u, 
		     LocalDataContainer<Scalar> const & v) {
      size_t n = u.getSize();
      if (n > v.getSize()) {
	RVLException e; e<<"Error: RVLLinComb::operator()\n";
	e<<"operand 1 longer than operand 2 - not enough data\n";
	throw e;
      }

      // compiler should generate more efficient code for 
      // several special cases by avoiding explicit multiply

      Scalar * pu = u.getData();
      Scalar const * pv = v.getData();

      // special case u = v + b*u
      if (a == ScalarFieldTraits<Scalar>::One()) {
	// special case u = v + u;
	if (b == ScalarFieldTraits<Scalar>::One()) {
	  for (size_t i=0;i<n;i++) {
	    pu[i] = pv[i]+pu[i];
	  }
	}
	// special case u = v - u;
	else if (b == -ScalarFieldTraits<Scalar>::One()) {
	  for (size_t i=0;i<n;i++) {
	    pu[i] = pv[i]-pu[i];
	  }
	}
	// special case u = v;
	else if (b == ScalarFieldTraits<Scalar>::Zero()) {
	  RVLCopy<Scalar> cp;    
	  // use RVLCopy which uses memcpy for efficiency
	  cp(u,v);
	}
      
	else { 
	  for (size_t i=0;i<n;i++) {
	    pu[i] = pv[i]+b*pu[i];
	  }
	}
      }
      // special case u = -v + b*w;
      else if (a == -ScalarFieldTraits<Scalar>::One()) {
	// special case u = -v;
	if (b == ScalarFieldTraits<Scalar>::Zero()) {
	  for (size_t i=0;i<n;i++) {
	    pu[i] = -pv[i];
	  }
	}
	else {
	  for (size_t i=0;i<n;i++) {
	    pu[i] = -pv[i]+b*pu[i];
	  }
	}
      }
      // special case u = a*v;
      else if (b == ScalarFieldTraits<Scalar>::Zero()) {
	for (size_t i=0;i<n;i++) {
	  pu[i] = a*pv[i];
	}
      }
      // general case u = a*v + b*u
      else {
	for (size_t i=0;i<n;i++) {
	  pu[i] = a*pv[i]+b*pu[i];
	}
      }
    }

    string getName() const { string s = "RVLLinCombObject"; return s; }
  };

  /** The standard linear algebra package which should work for most Spaces.
      Combines an elementwise assignment to zero, standard dot product, and elementwise
      linear combination.
  */
  template<class Scalar>
  class RVLLinearAlgebraPackage: public LocalLinearAlgebraPackage<Scalar,Scalar> {

  private:

    mutable RVLAssignConst<Scalar> this_zero;
    mutable RVLL2innerProd<Scalar> this_inner;
    mutable RVLLinCombObject<Scalar> this_lco;

  public:

    RVLLinearAlgebraPackage(Scalar ipscale = ScalarFieldTraits<Scalar>::One())
      : this_zero(ScalarFieldTraits<Scalar>::Zero()), 
	this_inner(abs(ipscale)), this_lco() {}
    RVLLinearAlgebraPackage(const RVLLinearAlgebraPackage<Scalar> & p) 
      :this_zero(ScalarFieldTraits<Scalar>::Zero()),
       this_inner(p.this_inner.getScale()),this_lco() {}
    ~RVLLinearAlgebraPackage() {}

    BinaryLocalFunctionObjectScalarRedn<Scalar, Scalar> & localinner() const {
      return this_inner; 
    }
    UnaryLocalFunctionObject<Scalar> & localzero() const {
      return this_zero;
    }
    LinCombObject<Scalar> & linComb() const {
      return this_lco; 
    }

    virtual bool compare(LinearAlgebraPackage<Scalar> const & lap) const {
      try {
	dynamic_cast<RVLLinearAlgebraPackage<Scalar> const &>(lap);
	return true;
      }
      catch (bad_cast) {
	return false;
      }
    }

    /** added to spparate instantiation from initialization */
    void setScale(Scalar newscale) { this_inner.setScale(newscale); }

    virtual void write(RVLException & e) const {
      e<<"RVLLinearAlgebraPackage\n";
    }

    virtual ostream & write(ostream & str) const {
      str<<"RVLLinearAlgebraPackage\n";
      return str;
    }
  };

}

#endif
