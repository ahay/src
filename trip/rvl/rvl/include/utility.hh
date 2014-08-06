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

#ifndef __RVL_UTILITY__
#define __RVL_UTILITY__

#include "except.hh"

/** Several of the main classes have operator new disabled. That is
    for good reason, see the TOMS paper, but some simply will not take
    no for an answer. If you must have operator new, define this
    macro.
 */
//#define RVL_OPERATOR_NEW_ENABLED

/** Include external random number generator to free unit testing from
    tyranny of compiler vendors. */
extern "C" {
#include "rngs.h"
}

/** The numeric precision is intended to convert type information into an 
    integer which we can use to control choices when interfacing with low--level code.
*/
template<class real>
inline static int numeric_precision() { return 0; }

template<>
inline int numeric_precision<float>() { return 1; }

template<>
inline int numeric_precision<double>() { return 2; }

namespace RVL {

  /**  A traits class to extend the capabilities of numeric_limits 
       without duplicating such capabilities.  This class and its
       specializations should have minimal methods and members

       ScalarType:  the redefined template type, as recommended by traits gurus
       AbsType:    the underlying real type, as would be returned by taking the absolute value
       Zero():        the additive identity
       One():         the multiplicative identity
     
       Note from ADP: I believe that computational methods (abs, conj,
       plus, minus, ect..) do not belong in this traits class.  They
       are well defined elsewhere and redirecting all computations
       through this class would only add unecessary inefficiency.

       I have borrowed some inspiration from the Teuchos::ScalarTraits struct,
       but 
       1. feel that struct is not the best approach 
       2. has a license we cannot use.
  */
  template<class Scalar>
  struct ScalarFieldTraits {
    typedef Scalar ScalarType;
    typedef Scalar AbsType;
    static inline Scalar Zero();
    static inline Scalar One();
    static inline Scalar AbsZero();
    static inline Scalar AbsOne();
  };

  template<>
  struct ScalarFieldTraits<bool> {
    typedef bool ScalarType;
    typedef bool AbsType;
    static inline int Zero() { return false; }
    static inline int One() { return true; }
    static inline int AbsZero() { return false; }
    static inline int AbsOne() { return true; }
  };

  template<>
  struct ScalarFieldTraits<int> {
    typedef int ScalarType;
    typedef int AbsType;
    static inline int Zero() { return 0; }
    static inline int One() { return 1; }
    static inline int AbsZero() { return 0; }
    static inline int AbsOne() { return 1; }
  };

  template<>
  struct ScalarFieldTraits<long> {
    typedef long ScalarType;
    typedef long AbsType;
    static inline long Zero() { return 0; }
    static inline long One() { return 1; }
    static inline long AbsZero() { return 0; }
    static inline long AbsOne() { return 1; }
  };

  template<>
  struct ScalarFieldTraits<unsigned int> {
    typedef unsigned int ScalarType;
    typedef unsigned int AbsType;
    static inline unsigned int Zero() { return 0; }
    static inline unsigned int One() { return 1; }
    static inline unsigned int AbsZero() { return 0; }
    static inline unsigned int AbsOne() { return 1; }
  };

  template<>
  struct ScalarFieldTraits<unsigned long> {
    typedef unsigned long ScalarType;
    typedef unsigned long AbsType;
    static inline unsigned long Zero() { return 0; }
    static inline unsigned long One() { return 1; }
    static inline unsigned long AbsZero() { return 0; }
    static inline unsigned long AbsOne() { return 1; }
  };


  template<>
  struct ScalarFieldTraits<float> {
    typedef float ScalarType;
    typedef float AbsType;
    static inline float Zero() { return 0.0; }
    static inline float One() { return 1.0; }
    static inline float AbsZero() { return 0.0; }
    static inline float AbsOne() { return 1.0; }
  };

  template<>
  struct ScalarFieldTraits<double> {
    typedef double ScalarType;
    typedef double AbsType;
    static inline double Zero() { return 0.0; }
    static inline double One() { return 1.0; }
    static inline double AbsZero() { return 0.0; }
    static inline double AbsOne() { return 1.0; }
  };

  template<class T>
  struct ScalarFieldTraits<std::complex<T> > {
    typedef complex<T> ScalarType;
    typedef T AbsType;
    static inline complex<T> Zero() { 
      return complex<T>(ScalarFieldTraits<T>::Zero());
    }
    static inline complex<T> One() {
      return complex<T>(ScalarFieldTraits<T>::One());
    }
    static inline T AbsZero() { return ScalarFieldTraits<T>::Zero(); }
    static inline T AbsOne() { return ScalarFieldTraits<T>::One(); }
  };

  /** RVL definition of complex conjugation. Uses std::conj for complex types,
      extends definition to basic floating point types so that uniform code
      can be written covering both real and complex cases.
  */
  static inline float conj(float x) { return x; }
  static inline double conj(double x) { return x; }
  static inline std::complex<float> conj(std::complex<float> x) { return std::conj(x); } 
  static inline std::complex<double> conj(std::complex<double> x) { return std::conj(x); } 

  /** Some applications do not make sense unless the absolute value
      type (signed, in the current implementation) is the same as the
      scalar type. Including a call to this function anywhere in the
      code defining a class will ensure that it compiles only when
      this "reality" condition is satisfied.
   */
  template<typename Scalar>
  void testRealOnly() {
    typename ScalarFieldTraits<Scalar>::AbsType a;
    Scalar b = ScalarFieldTraits<Scalar>::One();
    a=b;
    b=a;
  }
    
  /**  Calculate \f$quot = a/b\f$ in a careful manner.
       Without the tolerance, performs checks to avoid underflow/overflow.
       With the tolerance, checks to ensure that the resulting quotient 
       exceeds the tolerance.
       Return codes:
       1  Overflow
       2  Underflow
       3  Failed to exceed specified tolerance.
  */
  template<class real> 
  inline int ProtectedDivision(real a, real b, real & quot, 
			       real tol = ScalarFieldTraits<real>::AbsZero() ) {
    typedef typename ScalarFieldTraits<real>::AbsType absreal; 
    absreal fb = abs(b);
    absreal fa = abs(a);
    if (tol == ScalarFieldTraits<real>::AbsZero()) {
      if( fb < abs(ScalarFieldTraits<real>::One()) ) {
	if( fa < fb*std::numeric_limits<absreal>::max() ) {
	  quot = a/b;
	  return 0;
	} else {
	  return 1;
	}
      } else { // fb > 1.0
	if( fa > fb*numeric_limits<absreal>::min() ) {
	  quot = a/b;
	  return 0;
	} else {
	  quot = ScalarFieldTraits<real>::Zero();
	  //	return 2;
	  return 0;
	}
      }
    } else if( fb > fa*abs(tol) ) {
      quot = a/b;
      return 0;
    }
    else
      return 3;
  }

  /** Mixin interface to mandate write-to-ostream method, and derive
      write-to-exception method from it. Thanks to R. Bartlett. */
  class Writeable {

  public:
  
    /** Report state of object to ostream. */
    virtual ostream & write(ostream & str) const = 0;

      virtual ~Writeable() {};

    /** Report state of object to RVLException. */
    void write(RVLException & e) const {
      std::ostringstream ss;
      this->write(ss);
      e<<ss.str();
    }
  };

  /** Generic oracle interface. Should be able to say what it does, so
      child of Writeable. */
  template<typename T>
  class Oracle: public Writeable {
    
  public:

    /** yes or no */
    virtual bool isFeasible(T const &) const = 0;

  };

  /** Standard factory interface - really a policy. Should be able to
      say what it does, hence child of Writeable. */
  template<typename T>
  class Factory: public Writeable {

  public:
    
    Factory() {}
    Factory(Factory<T> const &) {}
    virtual ~Factory() {}

    /** generic build method */
    virtual T * build() const = 0;

  };


}

#endif
