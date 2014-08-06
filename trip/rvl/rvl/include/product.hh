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

#ifndef __RVL_PROD
#define __RVL_PROD

#include "utility.hh"
#include "except.hh"

namespace RVL {

  /** Simple Cartesian product mixin. */

  template<class T>
  class Product {

  public:
    
    Product() {}
    Product(const Product<T> &) {}
    virtual ~Product() {}

    virtual size_t getSize() const = 0;

    virtual T & operator[](size_t i) = 0;
    virtual T const & operator[](size_t i) const = 0;

  };

  /** Simple Cartesian product mixin - read-only access to factors. */

  template<class T>
  class ROProduct {

  public:
    
    ROProduct() {}
    ROProduct(const ROProduct<T> &) {}
    virtual ~ROProduct() {}

    virtual size_t getSize() const = 0;
    virtual T const & operator[](size_t i) const = 0;

  };
}

#endif
