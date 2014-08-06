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

#ifndef __RVL_LA
#define __RVL_LA

#include "data.hh"
#include "write.hh"

namespace RVL {

  /** LinCombObject is an subtype FunctionObject offering a method to
      set two mutable scalars a and b, intended to be coefficients in
      a linear combination. In reality the interface simply specifies
      accepting two scalars.*/

  template<class Scalar>
  class LinCombObject: public FunctionObject {
  public: 
    LinCombObject() {}
    LinCombObject(const LinCombObject<Scalar> &) {}
    virtual ~LinCombObject() {}
    /** Set method for linear combination coefficients */  
    virtual void setScalar(Scalar a, Scalar b) = 0;
  };

  /** Interface to function objects defining the basic ops
      of linear algebra in Hilbert space: linear combination,
      assignment to the zero vector, and inner product. Pure
      virtual base, permitting any convenient definition of
      these objects.      
  */

  template<class Scalar>
  class LinearAlgebraPackage: public Writeable {

  public:
    LinearAlgebraPackage() {}
    LinearAlgebraPackage(const LinearAlgebraPackage<Scalar> &) {}
    virtual ~LinearAlgebraPackage() {}
  
    /** access to inner product */
    virtual FunctionObjectScalarRedn<Scalar> & inner() const = 0;
    /** access to zero assignment */
    virtual FunctionObject & zero() const = 0;
    /** access to linear combination */
    virtual LinCombObject<Scalar> & linComb() const = 0;

    /** Compare for compatibility with another LinearAlgebraPackage.
	Usual comparison basis - is the type the same? However
	"compatibility" can be defined more loosely when appropriate:
	the intended meaning is "produces the same results when
	applied to the same data". Returns zero if not compatible,
	nonzero otherwise. 
    */
    virtual bool compare( LinearAlgebraPackage<Scalar> const & lap) const  = 0;

  };

}

#endif
